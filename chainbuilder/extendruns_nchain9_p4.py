# ingest the orbits*.h5 files to a single HDF5 file library/database
# 

import numpy as np
import sys
import os
import os.path
import h5py

import chainmodels

allmod =  {**chainmodels.runs_ext1 , **chainmodels.runs_ext2}
allfiles =  {**chainmodels.filenames_ext1 , **chainmodels.filenames_ext2}
for key in allmod :
  print(key,' ',allmod[key])

for key in chainmodels.runs :
  if chainmodels.runs[key][0] == 9 and chainmodels.runs[key][1] == 4:
    allmod[key] = chainmodels.runs[key]
    allfiles[key] = chainmodels.filenames[key]
    print(key,' ',allmod[key])

for key in allmod:
  print(key,' ',allmod[key])
  orb = h5py.File(allfiles[key], 'r+')

  #orb.attrs['tmax'] += 1e9*chainmodels.tkepler
  # this should work as the outputs arrys will just be allocated longer
  orb.attrs['nout'] = int((orb.attrs['tmax']+orb.attrs['dtout'])/orb.attrs['dtout']) +2
  orb.close()
