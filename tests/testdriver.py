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

driver = chaincalc.WallClockLimitedDriver(wall_start)

driver.runModel(tm, targettime, wall_limit, 100*keplertime)

tm2 = col.get_model(runi)
tm2.get_status()
print('Status from disc is ',tm2.status['status'])
if tm2.status['sim.t'] >= targettime:
  print('Run is evolved past targettime')


