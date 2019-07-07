# calls custom Rebound setup, as a library, which outputs to a HDF5 file

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

nchain = sys.argv[1]
p = sys.argv[2]
seqnum = int(sys.argv[3])

# run for 1e8 Kepler times
pert = 1e4*seqnum # vary the end time of ramping
keplertime = 2.0 * np.pi * 0.1**1.5 # orbit at 0.1
tdep = 2e5 * keplertime
deltatdep = tdep + pert
tmax = tdep + deltatdep + 1e8 * keplertime

# fire it off in pure C
run_sim(9, 6, tmax, tdep, deltatdep, 125)
