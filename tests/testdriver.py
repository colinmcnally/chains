# Colin McNally 2019 <colin@colinmcnally.ca>
# Script to test full usage.
# Set to run in laptop, not container
# 
import sys
sys.path.append('/Users/colinm/vc/rebound')
sys.path.append('../crbx')
sys.path.append('../')
import time
import rebound
import chaincalc

wall_start = time.time()

wall_limit = 25

from testcampaign import *

#Now have name col from testcampaign

runi = int(sys.argv[1])
print('Will try to run model at index {} of testcampaign'.format(runi))

tm = col.get_model(runi)

print('hash is ',tm.hash)

tm.init_rebound()

print('status after init')
tm.sim.status();

if tm.status['status']=='running':
    tm.set_wall_start(wall_start)

    try:
        lastsnap = wall_start
        while (tm.sim.t < targettime) and ((time.time() - wall_start) <= wall_limit) and tm.status['status']=='running':
            nextcheck = min(tm.sim.t+100*keplertime, targettime)
            tm.sim.integrate(nextcheck)
            if (time.time() - lastsnap) > wall_check_interval:
                print('\nmanual snapshot wall elapsed:', time.time()-wall_start)
                tm.manual_snapshot()
                lastsnap = time.time()
        if (time.time() - wall_start) > wall_limit:
            print('\nwall limit exceeded:', time.time()-wall_start)
            tm.manual_snapshot()
    except rebound.Collision:
        tm.mark_collision()
    else:
        tm.manual_snapshot()
    print('')
else:
    print('Not integrating, the status is already {}'.tm.status['status'])


tm2 = col.get_model(runi)
tm2.get_status()
print('Status from disc is ',tm2.status['status'])
if tm2.status['sim.t'] >= targettime:
  print('Run is evolved past targettime')


