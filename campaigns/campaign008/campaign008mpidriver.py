# Colin McNally 2019 <colin@colinmcnally.ca>
# MPI based driver for campaign008 
#
import time
import chaincalc
import rebound

print("Testing version - bound from disc")

wall_start = time.time()

# three days minus a bit on apocrita
wall_limit_total = 72*60*60 - 15*60

#Intervals of wallclock advance of a given model
wall_limit_chunk = 27*60

from campaign008 import *
col = Campaign008()
# now have name col 

check_interval = 100.0*keplertime

class CheckChainSheduler(chaincalc.MpiSchedulerAsync):
    """set the targettime_hook method"""

    def targettime_hook(self, model):
        #This gets called whenever a model is being run, even if not 
        # evolved because it was already finished, so avoid doing this twice.
        findingdkey = 'is_this_a_constant_ratio_chain'
        if not (findingdkey in model.status):
            sa = rebound.SimulationArchive(model.status['simarchive_filename'], process_warnings=False)
            oa = chaincalc.OrbitArray(sa)
            tlook = 1e4*keplertime #probably 1e2-1e3 times the output cadence
            lko = oa.tail_orbitarray(tlook)
            finding = oa.is_this_a_constant_ratio_chain(oa.compute_tight_angles(lko))
            model.status['is_this_a_constant_ratio_chain'] = finding
            model.overwrite_status()


#init does the running...
runner = CheckChainSheduler(col, wall_start, wall_limit_total, check_interval, wall_limit_chunk=wall_limit_chunk)
