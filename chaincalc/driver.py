# Colin McNally <colinmcnally.ca>
# Base classes for drivers
# The first implemented will be an abstraction of the existing job array-based driver scripts
# Should also design a MPI-parallel master-slave system
#
import time
import rebound

class DriverBase():
    def __init__(self, physical_outputs):
        self.physical_outputs = sorted(physical_outputs) 
  
    def runModel(self, model):
        pass
 
    def set_physical_outputs(self, physical_outputs):
       self.physical_outputs = physical_outputs

    def get_next_physical_output(self, t):
       """Return the next physical output time past time t, assumes a list of times, not a frequency, might be a limitation"""
       if (len(self.physical_outputs) > 0):
           #warning, this will not scale too far, need to implement a list frequency option
           for pt in self.physical_outputs:
               if pt > t:
                   return pt
           return float('inf')
       else:
           return 0.0

class WallClockLimitedDriver(DriverBase):
    """An abstracted driver class for wallclock-limited execution."""

    def __init__(self, wall_start):
        """Initialize driver, and keep track of the process start time wall_start, calls super class __init__"""
        self.wall_start = wall_start


    def runModel(self, model, targettime, wall_limit, check_interval):
        """Run the model passed, wall_start should be the job start time, and check execution on physical time check_time"""

        self.set_physical_outputs(model.params['physical_outputs'])
        model.init_rebound()

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
                        print('\nManual snapshot flag: {} wall elapsed: {:e}'.format(manual_flag, time.time() -self.wall_start))
                        model.manual_snapshot()
                        lastsnap = time.time()

                if (time.time() - self.wall_start) > wall_limit:
                    print('\nWall limit exceeded, will checkpoint and shutdown: {:e}'.format(time.time() -self.wall_start))
                    model.manual_snapshot()

            except rebound.Collision:
                model.mark_collision()
            else:
                model.manual_snapshot()
            model.unlock()
            print('')
        else:
            print('Not integrating, the status is already: {}',model.status['status'])
