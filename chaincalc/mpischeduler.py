# Colin McNally 2019 <colin@colinmcnally.ca>
# A class that orchastrates campaign advances in fixed wallclock chunks
# 
from mpi4py import MPI
import collections
import time
import chaincalc

class MpiScheduler:
    """A basic MPI based schedualer that works in wallclock chunks. All processes do this together in a synchronous manner."""
    def __init__(self, campaign, wall_start, wall_limit_total, check_interval, wall_limit_chunk = 5*60.0):
        """__init__ calls everything, this is execute-on-instantiate"""
        self.exitFlag = False
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nranks = self.comm.Get_size()
        self.wall_start = wall_start
        self.wall_limit_total = wall_limit_total
        self.wall_limit_chunk = wall_limit_chunk
        if self.rank == 0:
             self.masterprocess(campaign, check_interval)
        else:
             self.slaveprocess(campaign, check_interval)


    def masterprocess(self, campaign, check_interval):
        """This is the master process, make a work list and send out assignments, do some work too"""

        print('Master Process rank ',self.rank)
        am = collections.deque(range(0, campaign.get_size()))
        finishedruns = []
        # go into a loop sending out and doing work, use blocking comms so I/O is regulated a bit
        exitFlag = False
        while (time.time() < self.wall_start + self.wall_limit_total - self.wall_limit_chunk and
               not exitFlag):

            # send work to processes, only send out as much work as we have
            for dest in range(1,self.nranks):
                if len(am) > 0: 
                    torun = am.popleft()
                    self.comm.send(['run', torun], dest=dest, tag=1)
                else:
                    self.comm.send(['wait'], dest=dest, tag=1)

            # if there is work leftover, take it
            if len(am) > 0: 
                runi = am.popleft()
                rundone = self.run_model(campaign, runi, check_interval)
                # could make this return more symmetrical
                if rundone:
                    finishedruns.append(runi)
                else:
                    am.append(runi)

            # sent soemthing to all slaves, should get something back
            for src in range(1,self.nranks):
                finished = self.comm.recv(source=src, tag=2) 
                if finished[0] == 'finished':
                    finishedruns.append(finished[1])
                elif finished[0] == 'running':
                    #add this run to back of queue
                    am.append(finished[1])
                elif finished[0] == 'waiting':
                    pass
                else:
                    print('Error, MpiScheduler unknown status message ',finished)
            print('***Master working list is now ',list(am))
            if len(am) == 0:
                for dest in range(1,self.nranks):
                    self.comm.send(['exit'], dest=dest, tag=1)
                exitFlag = True
        print('Exit rank ',self.rank)


    def slaveprocess(self, campaign, check_interval):
        """Execution of a MPI slave process, waits for commands from master, executes chunks of work from campaign"""
        print('Slave Process rank ',self.rank)

        exitFlag = False
        while (time.time() < self.wall_start + self.wall_limit_total - self.wall_limit_chunk and not exitFlag):
            #post a blocking recv, to wait for an assignment from rank==0
            cmd = self.comm.recv(source=0, tag=1)
            if cmd[0] == 'run': 
                rundone = self.run_model(campaign, cmd[1], check_interval)
                if rundone:
                    self.comm.send(['finished',cmd[1]], dest=0, tag=2)
                else:
                    self.comm.send(['running',cmd[1]], dest=0, tag=2)
            elif cmd[0] == 'exit':
                exitFlag = True
            else:
                self.comm.send(['waiting'], dest=0, tag=2)
        print('Exit rank ',self.rank)


    def run_model(self, campaign, runi, check_interval):
        """Execute a model from the campaign, return True if it has passed targettime"""

        tm = campaign.get_model(runi)
        targettime = campaign.get_targettime(runi)
        print('Rank {} will try to run model at index {} of campaign hash {}'.format(self.rank, runi, tm.hash), flush=True)
        now = time.time()
        driver = chaincalc.WallClockLimitedDriver(now)
        # run a chunk, but don't run past the total lifetime allowed for this process
        now_limit = min(self.wall_start +self.wall_limit_total, now +self.wall_limit_chunk) - now
        driver.runModel(tm, targettime, now_limit, check_interval)

        tm2 = campaign.get_model(runi)
        tm2.get_status()
        if tm2.status['sim.t'] >= targettime:
            print('Run is evolved past targettime')
            self.targettime_hook(tm2)
            return True
        elif tm2.status['lock']:
            return True
        else:
            return False

    def targettime_hook(self, model):
        """A method to be called once a model has exceeded targettime"""
        pass 


class MpiSchedulerAsync(MpiScheduler):
    """ This version of the MPI schedualer lets the root process only manage work, and slaves get work whenever they need it."""

    def masterprocess(self, campaign, check_interval):
        """This is the master process, make a work list and send out assignments, do not work."""

        print('Master Process rank ',self.rank, flush=True)
        am = collections.deque(range(0, campaign.get_size()))
        finishedruns = []
        waiting = []
        
        # loop over workers, assign everything a chunk of work if possible, and if not assign wait
        # This can be blocking sends.
        # send work to processes, only send out as much work as we have
        if (time.time() < self.wall_start + self.wall_limit_total - self.wall_limit_chunk):
            for dest in range(1,self.nranks):
                if len(am) > 0: 
                    torun = am.popleft()
                    self.comm.send(['run', torun], dest=dest, tag=1)
                else:
                    self.comm.send(['wait'], dest=dest, tag=1)

        # go into a loop sending out and doing work
        exitFlag = False
        while (time.time() < self.wall_start + self.wall_limit_total - self.wall_limit_chunk and
               not exitFlag):
        
            # Event handling loop until close to wallclock limit
            # set a blocking probe for incoming with tag=2
            # get a return message
            # if model finished, put it on the finished list
            # else if model not finished, add it back to queue
            # if the process was waiting, let it keep waiting, it will be in a blocking recv
            # if no more work to do, send exit command to all processes, and set out exitFlag.

            status = MPI.Status()
            s = self.comm.probe(source=MPI.ANY_SOURCE, tag=2, status=status)
            #print('Root Probe got ', s, status.source, flush=True)
            finished = self.comm.recv(source=status.source, tag=2) 
            if finished[0] == 'finished':
                finishedruns.append(finished[1])
            elif finished[0] == 'running':
                #add this run to back of queue
                am.append(finished[1])
            elif finished[0] == 'waiting':
                if status.source not in waiting:
                    waiting.append(status.source)
            else:
                print('Error, MpiScheduler unknown status message ',finished)
            
            #reply address
            dest = status.source
            # send work to process, or tell it to wait (only once!)
            if len(am) > 0: 
                torun = am.popleft()
                self.comm.send(['run', torun], dest=dest, tag=1)
            elif dest not in waiting:
                # important not to saed wait twice, a waiting dest rank just needs an exit
                self.comm.send(['wait'], dest=dest, tag=1)

            print('***Master unassigned list is now ',list(am), flush=True)
            if len(am) == 0 and len(waiting) == self.nranks-1:
                exitFlag = True

        #in this Async mode we need to make sure to send this to slaves waiting on the recv before exiting
        for dest in range(1,self.nranks):
            # we never track this, just let it die
            self.comm.isend(['exit'], dest=dest, tag=1)
        print('Exit rank ',self.rank, flush=True)


    def slaveprocess(self, campaign, check_interval):
        """Execution of a MPI slave process, waits for commands from master, executes chunks of work from campaign"""
        print('Slave Process rank ',self.rank, flush=True)

        exitFlag = False
        while (time.time() < self.wall_start + self.wall_limit_total - self.wall_limit_chunk and not exitFlag):
            #post a blocking recv, to wait for an assignment or exit from rank==0
            cmd = self.comm.recv(source=0, tag=1)
            if cmd[0] == 'run': 
                rundone = self.run_model(campaign, cmd[1], check_interval)
                if rundone:
                    # These are all non-blocking so we can get interrupted by a shutdown if needed
                    request = self.comm.isend(['finished',cmd[1]], dest=0, tag=2)
                else:
                    request = self.comm.isend(['running',cmd[1]], dest=0, tag=2)
            elif cmd[0] == 'exit':
                exitFlag = True
            else:
                request = self.comm.isend(['waiting'], dest=0, tag=2)
        print('Exit rank ',self.rank, flush=True)


