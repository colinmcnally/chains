# define a model grid
# nchain - from 3 to 10
# p - the resonance in p+1:p 4:3 to 7:6
# realization - say 10 realizations of each
import numpy as np
from hashlib import sha1

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

realizationaxis_ext1 = np.arange(10, 50, dtype=np.int)
p = 4
ncahin = 4
runs_ext1 = {}
filenames_ext1 = {}
keys_ext1 = {}
for realization in realizationaxis:
  runtuple = (nchain, p, realization)
  ha = sha1()
  ha.update(bytes(runtuple))
  key = ha.digest().hex()
  runs_ext1[key] = runtuple
  keys_ext1[runtuple] = key
  filenames_ext1[key] = 'orbits_{}_{}:{}_{:.2e}_{}.h5'.format(nchain, p+1, p, 1e-5, realization) 
 
def merge_dict(a,b):
  c = {}
  for d in (a, c):
    c.update(d)
  return c 

runs_all = merge_dict(runs, runs_ext1)
filenames_all =  merge_dict(filenames, filenames_ext1)
keys_all =  merge_dict(keys, keys_ext1)
