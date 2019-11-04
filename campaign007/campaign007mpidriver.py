# Colin McNally 2019 <colin@colinmcnally.ca>
# MPI based driver for campaign007 
#
import time
import chaincalc

wall_start = time.time()

# three days minus a bit on apocrita
wall_limit_total = 72*60*60 - 15*60

#Intervals of wallclock advance of a given model
wall_limit_chunk = 27*60

from campaign007 import *
col = Campaign007()
# now have name col 

check_interval = 100.0*keplertime

#init does the running...
runner = chaincalc.MpiScheduler(col, targettime, wall_start, wall_limit_total, check_interval, wall_limit_chunk=wall_limit_chunk)
