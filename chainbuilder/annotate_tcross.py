# build a Pandas hdf5 of the system lifetime table for faster usage in plotting

import chainmodels
import numpy as np
import numpy.lib.recfunctions as rfn
import sys
import h5py
import matplotlib.pyplot as plt
import scipy.stats
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.neighbors.kde import KernelDensity
from sklearn.covariance import EllipticEnvelope, MinCovDet
import pandas as pd
import scipy.integrate


tkepler = chainmodels.tkepler

def build_lifetable(rundict, filename):
  '''
  Builds tables of system lifetimes, etc from a library file
  '''
  library_file = h5py.File( filename, 'r')
  lifetimes = {}
  lasttimes = {}
  for key in rundict:
    orb = library_file['/'+key]
    tcoll = orb.attrs['tcollstop'][0]
    tstart = orb.attrs['tdep'][0]+orb.attrs['deltatdep'][0]
    if tcoll > 0.0:
      tlife = tcoll - tstart
    else:
      tlife = -1 
    lifetimes[key] = tlife/tkepler
    lasttimes[key]  = (orb['t'].value.max()-tstart)/tkepler

  library_file.close()

  lifetable = np.zeros(len(rundict), 
                     dtype=[('nchain',np.int16), ('p',np.int16),
                            ('realization',np.int16), ('tlife',np.float64),
                            ('lasttime',np.float64)])
  for i,key in enumerate(rundict):
    nchain, p, realization = rundict[key]
    lifetable[i]['nchain'] = nchain
    lifetable[i]['p'] = p
    lifetable[i]['realization'] = realization
    lifetable[i]['tlife'] = lifetimes[key]
    lifetable[i]['lasttime'] = lasttimes[key]

  models = set()
  for i,key in enumerate(rundict):
    nchain, p, realization = rundict[key]
    models.add((nchain,p))

  return [models, lifetable]

# build tables of the system lifetimes form library files
filename = 'chain_data_library3.h5'
models_ext2, lifetable_ext2 = build_lifetable(chainmodels.runs_ext2, filename)
filename = 'chain_data_library2.h5'
models_ext1, lifetable_ext1 = build_lifetable(chainmodels.runs_ext1, filename)
filename = 'chain_data_library.h5'
models_old, lifetable_old = build_lifetable(chainmodels.runs, filename)


lifetable = rfn.stack_arrays((lifetable_ext2, lifetable_ext1, lifetable_old), usemask=False)
print('merged lifetable shape ',lifetable.shape, lifetable.dtype)
models = set().union(models_ext2).union(models_ext1).union(models_old)

lifetableframe = pd.DataFrame(lifetable)
modelsframe = pd.DataFrame(list(models), columns=['nchain','p'])
lifetableframe.to_hdf('lifetable.h5','lifetable')
modelsframe.to_hdf('lifetable.h5','models')

