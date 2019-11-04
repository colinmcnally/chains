# Colin McNally 2019 <colin@colinmcnally.ca>
# A class that orchatrates campaign  advances in fixed wallclock chunks, using round-robin schedualing
# 
from mpi4py import MPI
import collections
import time
import chaincalc

class MpiScheduler:
    def __init__(self, campaign, targettime, wall_start, wall_limit_total, check_interval, wall_limit_chunk = 5*60.0):
        """__init__ calls everything, this is execute-on-instantiate"""
        self.exitFlag = False
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nranks = self.comm.Get_size()
        self.wall_start = wall_start
        self.wall_limit_total = wall_limit_total
        self.wall_limit_chunk = wall_limit_chunk
        if self.rank == 0:
             self.masterprocess(campaign, targettime, check_interval)
        else:
             self.slaveprocess(campaign, targettime, check_interval)


    def masterprocess(self, campaign, targettime, check_interval):
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
                rundone = self.run_model(campaign, runi, targettime, check_interval)
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
                    print('Error, MpiSchedualer unknown status message ',finished)
            print('***Master working list is now ',list(am))
            if len(am) == 0:
                for dest in range(1,self.nranks):
                    self.comm.send(['exit'], dest=dest, tag=1)
                exitFlag = True
        print('Exit rank ',self.rank)


    def slaveprocess(self, campaign, targettime, check_interval):
        """Execution of a MPI slave process, waits for commands from master, executes chunks of work from campaign"""
        print('Slave Process rank ',self.rank)

        exitFlag = False
        while (time.time() < self.wall_start + self.wall_limit_total - self.wall_limit_chunk and not exitFlag):
            #post a blocking recv, to wait for an assignment from rank==0
            cmd = self.comm.recv(source=0, tag=1)
            if cmd[0] == 'run': 
                rundone = self.run_model(campaign, cmd[1], targettime, check_interval)
                if rundone:
                    self.comm.send(['finished',cmd[1]], dest=0, tag=2)
                else:
                    self.comm.send(['running',cmd[1]], dest=0, tag=2)
            elif cmd[0] == 'exit':
                exitFlag = True
            else:
                self.comm.send(['waiting'], dest=0, tag=2)
        print('Exit rank ',self.rank)


    def run_model(self, campaign, runi, targettime, check_interval):
        """Execute a model from the campaign, return True if it has passed targettime"""

        tm = campaign.get_model(runi)
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
            return True
        else:
            return False


