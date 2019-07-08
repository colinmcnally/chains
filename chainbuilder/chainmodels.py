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
for nchain in nchainaxis:
  for p in paxis:
    for realization in realizationaxis:
      runtuple = (nchain, p, realization)
      ha = sha1()
      ha.update(bytes(runtuple))
      key = ha.digest().hex()
      runs[key] = runtuple
      filenames[key] = 'orbits_{}_{}:{}_{:.2e}_.h5'.format(nchain, p+1, p, 1e-5, realization) 
 
#output - sequence number - the hash of the runtuple

#needs to give a two-way mapping
