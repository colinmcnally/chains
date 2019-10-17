# ingest the orbits*.h5 files to a single HDF5 file library/database
# 

import numpy as np
import sys
import os
import os.path
import h5py

import chainmodels


libf = h5py.File('chain_data_library2.h5', 'w')

#sha1 gives 40 byte hash. so this datatype relies on that choice as the run encoding
runtable = np.zeros(len(chainmodels.runs_ext1), 
                    dtype = [('key', 'S40'), ('nchain', np.int), ('p', np.int), ('realization', np.int)])

# need to do a bit of gynmnastics around the coding/deconding to/from HDF5 through numpy
# to use the key reda back from h5pyi as a bytes object, use the bytes.hex() method to get a string matching the 
# one used as a group name below
runtable[:]['key'] = [ bytes.fromhex(key) for key in chainmodels.runs_ext1.keys()]

runtuples = np.asarray(list(chainmodels.runs_ext1.values()))
runtable[:]['nchain']      = runtuples[:,0]
runtable[:]['p']           = runtuples[:,1]
runtable[:]['realization'] = runtuples[:,2]

runsh5 = libf.create_dataset(b'runtable', data = runtable)
runsh5.attrs.create(b'headers', (b'nchain', b'p', b'realization'))

for key in chainmodels.runs_ext1:
  print(key,' ',chainmodels.runs_ext1[key])
  orb = h5py.File(chainmodels.filenames_ext1[key], 'r')

  nplanets = orb.attrs['nchain'][0]
  ie     = orb.attrs['lastout'][0]
  tcoll  = orb.attrs['tcollstop']
  tstart = orb.attrs['tdep'][0] +orb.attrs['deltatdep'][0]

  dt = orb['t'][ie-1] / (ie-1)
  mask = np.ones(orb['t'][:ie].shape, dtype=np.bool)
  times = orb['t'][0:ie]
  for i,t in enumerate(times[1:]):
    if (t<times[i]):
      mask[i+1] = False

  obg = libf.create_group(key)
  for att in orb.attrs.items():
    obg.attrs.create(att[0], att[1]) 

  for okey in orb.keys():
    if len(orb[okey].shape) >1:
      obg.create_dataset( okey, data = np.compress(mask, orb[okey][:,:ie], axis=1) )
    else:
      obg.create_dataset( okey, data = np.compress(mask, orb[okey][:ie], axis=0) )

  orb.close()
libf.close()
