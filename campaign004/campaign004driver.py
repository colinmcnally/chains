# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign004 integration driver
# 
import sys
import time
import rebound
import chaincalc

wall_start = time.time()


# three days minus a bit on apocrita
wall_limit = 72*60*60 - 15*60

from campaign004 import *

#Now have name col from campaign001
col = Campaign004()
print('Campaign has {} models'.format(col.get_size()))

# Job array indexes start at 1 not 0
runi = int(sys.argv[1])-1
print('Will try to run model at index {} of campaign004'.format(runi))

print('Wallclock checkpoint interval is {} seconds'.format(wall_check_interval))


tm = col.get_model(runi)

print('hash is ',tm.hash)

driver = chaincalc.WallClockLimitedDriver(wall_start)

driver.runModel(tm, targettime, wall_limit, 100*keplertime)

tm2 = col.get_model(runi)
tm2.get_status()
print('On exit, status from disc is ',tm2.status['status'])
if tm2.status['sim.t'] >= targettime:
  print('Run is evolved past targettime')


