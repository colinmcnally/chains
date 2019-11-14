# Colin McNally <colinmcnally.ca>
# Base classes for drivers
# The first implemented will be an abstraction of the existing job array-based driver scripts
# Should also design a MPI-parallel master-slave system
#
import time
import rebound

class DriverBase():
    def __init__(self, verbose=False):
        self.verbose = verbose
  
    def runModel(self, model):
        pass
 
    def set_physical_outputs(self, physical_outputs):
       self.physical_outputs = sorted(physical_outputs)

    def get_next_physical_output(self, t):
       """Return the next physical output time past time t, assumes a list of times, not a frequency, might be a limitation"""
       nt = float('inf')
       if (self.model.params['physical_otuput_dt'] is not None):
           nt = self.model.status.last_physical_output + self.model.params['physical_otuput_dt']
       if (len(self.physical_outputs) > 0):
           #warning, this will not scale too far, need to implement a list frequency option
           for pt in self.physical_outputs:
               if pt > t:
                   return min(nt, pt)
           return min(nt,float('inf'))
       else:
           return float('inf')

class WallClockLimitedDriver(DriverBase):
    """An abstracted driver class for wallclock-limited execution."""

    def __init__(self, wall_start, verbose=False):
        """Initialize driver, and keep track of the process start time wall_start, calls super class __init__"""
        self.wall_start = wall_start
        self.verbose = verbose


    def runModel(self, model, targettime, wall_limit, check_interval):
        """Run the model passed, wall_start should be the job start time, and check execution on physical time check_time"""

        # set model verbose following driver verbose
        model.verbose = self.verbose
        model.init_rebound()
        self.model = model

        self.set_physical_outputs(model.params['physical_outputs'],
                                  model.params['physical_output_dt'],
                                  model.status['last_physical_output'])

        if self.verbose:
            print('Model status after init')
            model.sim.status()

        if model.status['status']=='running':
            if model.status['lock']:
                print('Error, the model is locked. Will not modify see:',model.status_filename)
            else:
                model.lock()
                model.set_wall_start(self.wall_start)

                try:
                    lastsnap = self.wall_start

                    while (model.sim.t < targettime) and ((time.time() - self.wall_start) <= wall_limit) and model.status['status']=='running':
                        nextcheck = min(model.sim.t+check_interval, targettime)
                        nextphysout = self.get_next_physical_output(model.sim.t)
                        if nextcheck > nextphysout:
                            manual_flag = True    
                        else:
                            manual_flag = False   
         
                        model.sim.integrate(nextcheck)
                        if ((time.time() - lastsnap) > model.params['snap_wall_interval']) or manual_flag:
                            if self.verbose:
                                print('\nManual snapshot flag: {} wall elapsed: {:e}'.format(manual_flag, time.time() -self.wall_start))
                            model.manual_snapshot()
                            lastsnap = time.time()

                    if (time.time() - self.wall_start) > wall_limit:
                        if self.verbose:
                            print('\nWall limit exceeded, will checkpoint and shutdown: {:e}'.format(time.time() -self.wall_start))
                        model.manual_snapshot()

                except rebound.Collision:
                    model.mark_collision()
                else:
                    model.manual_snapshot()
                model.unlock()
                if self.verbose:
                    print('')
        else:
            if self.verbose:
                print('Not integrating, the status is already: {}',model.status['status'])
