# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign002 integration driver
# 
import sys
import time
import rebound
import chaincalc

wall_start = time.time()


# three days minus a bit on apocrita
wall_limit = 72*60*60 - 15*60

from campaign002 import *

col = Campaign002()
print('Campaign has {} models'.format(col.get_size()))

###This thing could all be put in the chaincalc module

# Job array indexes start at 1 not 0
runi = int(sys.argv[1])-1
print('Will try to run model at index {} of campaign002'.format(runi))



tm = col.get_model(runi)
print('Wallclock checkpoint interval is {} seconds'.format(tm.params['snap_wall_interval']))

print('hash is ',tm.hash)

tm.init_rebound()

print('status after init')
tm.sim.status();

if tm.status['status']=='running':
    if tm.status['lock']:
        print('Error, the model is locked. Will not modify see:',tm.status_filename)
    else:
        tm.lock()
        tm.set_wall_start(wall_start)

        try:
            lastsnap = wall_start
            while (tm.sim.t < targettime) and ((time.time() - wall_start) <= wall_limit) and tm.status['status']=='running':
                nextcheck = min(tm.sim.t+100*keplertime, targettime)
                tm.sim.integrate(nextcheck)
                if (time.time() - lastsnap) > tm.params['snap_wall_interval']:
                    print('\nmanual snapshot wall elapsed:', time.time()-wall_start)
                    tm.manual_snapshot()
                    lastsnap = time.time()
            if (time.time() - wall_start) > wall_limit:
                print('\nwall limit exceeded, will checkpoint and shutdown:', time.time()-wall_start)
                tm.manual_snapshot()
        except rebound.Collision:
            tm.mark_collision()
        else:
            tm.manual_snapshot()
        tm.unlock()
        print('')
else:
    print('Not integrating, the status is already {}'.tm.status['status'])


tm2 = col.get_model(runi)
tm2.get_status()
print('On exit, status from disc is ',tm2.status['status'])
if tm2.status['sim.t'] >= targettime:
  print('Run is evolved past targettime')


