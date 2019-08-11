# define a model grid
# nchain - from 3 to 10
# p - the resonance in p+1:p 4:3 to 7:6
# realization - say 10 realizations of each
import numpy as np
from hashlib import sha1

tkepler = 2*np.pi*0.1**1.5

nchainaxis = np.arange(3, 11, dtype=np.int)
paxis = np.arange(3, 7, dtype=np.int)
realizationaxis = np.arange(0, 10, dtype=np.int)

runs = {}
filenames = {}
keys = {}
for nchain in nchainaxis:
  for p in paxis:
    for realization in realizationaxis:
      runtuple = (nchain, p, realization)
      ha = sha1()
      ha.update(bytes(runtuple))
      key = ha.digest().hex()
      runs[key] = runtuple
      keys[runtuple] = key
      filenames[key] = 'orbits_{}_{}:{}_{:.2e}_{}.h5'.format(nchain, p+1, p, 1e-5, realization) 
 
#output - sequence number - the hash of the runtuple


#add extensions to runs for a single chain configuration
realizationaxis_ext1 = np.arange(10, 50, dtype=np.int)
p = 4
nchain = 9
runs_ext1 = {}
filenames_ext1 = {}
keys_ext1 = {}
for realization in realizationaxis_ext1:
  runtuple = (nchain, p, realization)
  ha = sha1()
  ha.update(bytes(runtuple))
  key = ha.digest().hex()
  runs_ext1[key] = runtuple
  keys_ext1[runtuple] = key
  filenames_ext1[key] = 'orbits_{}_{}:{}_{:.2e}_{}.h5'.format(nchain, p+1, p, 1e-5, realization) 


# this manner of adding extensions is sort of awkward, but it does allow for lots of generality
realizationaxis_ext2 = np.arange(50, 250, dtype=np.int)
p = 4
nchain = 9
runs_ext2 = {}
filenames_ext2 = {}
keys_ext2 = {}
for realization in realizationaxis_ext2:
  runtuple = (nchain, p, realization)
  ha = sha1()
  ha.update(bytes(runtuple))
  key = ha.digest().hex()
  runs_ext2[key] = runtuple
  keys_ext2[runtuple] = key
  filenames_ext2[key] = 'orbits_{}_{}:{}_{:.2e}_{}.h5'.format(nchain, p+1, p, 1e-5, realization) 
  
runs_all = {**runs, **runs_ext1, **runs_ext2}
filenames_all =  {**filenames, **filenames_ext1, **filenames_ext2}
keys_all =  {**keys, **keys_ext1, **keys_ext2}
