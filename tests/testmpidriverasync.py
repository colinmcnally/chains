# Colin McNally 2019 <colin@colinmcnally.ca>
# Demo of a MPI based driver enabling asyncronous evolution
# advances sims from a queue in wallclock chunks
# Set to run in laptop, not container
# 
#from mpi4py import MPI
import collections
import sys
sys.path.append('/Users/colinm/vc/rebound')
sys.path.append('../crbx')
sys.path.append('../')
import time
import rebound
import chaincalc


wall_start = time.time()

wall_limit_total = 60*5

#
wall_limit_chunk = 10


from testcampaign import *
# now have name col 

check_interval = 100.0*keplertime

#init does the running...
runner = chaincalc.MpiSchedulerAsync(col, targettime, wall_start, wall_limit_total, check_interval, wall_limit_chunk=wall_limit_chunk)
