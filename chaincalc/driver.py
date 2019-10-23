# Colin McNally <colinmcnally.ca>
# Base classes for drivers
# The first implemented will be an abstraction of the existing job array-based driver scripts
# Should also design a MPI-parallel master-slave system
#
import time
import rebound
#from . import chainmodel

class DriverBase():
    def runModel(self, model):
        pass


class WallClockLimitedDriver(DriverBase):
    """An abstracted driver class for wallclock-limited execution."""

    def __init__(self, wall_start):
        """Initialize driver, and keep track of the process start time wall_start"""
        self.wall_start = wall_start


    def runModel(self, model, targettime, wall_limit, check_interval):
        """Run the model passed, wall_start should be the job start time, and check execution on physical time check_time"""

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
                    model.sim.integrate(nextcheck)
                    if (time.time() - lastsnap) > model.params['snap_wall_interval']:
                        print('\nManual snapshot wall elapsed: {:e}'.format(time.time() -self.wall_start))
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
