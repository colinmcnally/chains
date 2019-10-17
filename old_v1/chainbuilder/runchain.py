# calls custom Rebound setup, as a library, which outputs to a HDF5 file

import chainmodels
import ctypes
import numpy as np
import sys


# load the chainbuilder library, and set up arrgtypes for the run_sim symbol
if sys.platform == 'linux': # assume our Singularity container
  _hdf5 = ctypes.CDLL('/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so', mode=ctypes.RTLD_GLOBAL)
  _hdf5_hl = ctypes.CDLL('/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.so', mode=ctypes.RTLD_GLOBAL)
  _reb = ctypes.CDLL('/opt/rebound/src/librebound.so', mode=ctypes.RTLD_GLOBAL)
  _chain = ctypes.CDLL('/opt/chains/chainbuilder/chainbuilder.so', mode=ctypes.RTLD_GLOBAL)
else: # probably macOS on my laptop
  _chain = ctypes.CDLL('chainbuilder.so')
_chain.run_sim.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int)

# a little wrapper
def run_sim( nchain, p, tmax, tdep, deltatdep, seqnum):
    global _chain
    _chain.run_sim(ctypes.c_int(nchain), ctypes.c_int(p), ctypes.c_double(tmax), ctypes.c_double(tdep), 
                   ctypes.c_double(deltatdep), ctypes.c_int(seqnum))

#set run times
runindex = int(sys.argv[1])-1 #job array indexes start at 1

keys = list(chainmodels.runs.keys())
keys = sorted(keys)
hashkey = keys[runindex]
nchain, p, seqnum = chainmodels.runs[hashkey] 
print('chainmodel dict length ', len(chainmodels.runs))
print('nchain ', nchain,' p ',p,' realiz ', seqnum,' hashkey ',hashkey)

# run for 1e8 Kepler times
keplertime = 2.0 * np.pi * 0.1**1.5 # orbit at 0.1
pert = 1e3*keplertime*seqnum # vary the end time of ramping
tdep = 4e5 * keplertime + pert
deltatdep = tdep 
tmax = tdep + deltatdep + 1e8 * keplertime
if (nchain==9) and (p==4):
  tmax = tdep + deltatdep + 1e9 * keplertime


# fire it off in pure C
run_sim(nchain, p, tmax, tdep, deltatdep, seqnum)
