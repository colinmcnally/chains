# This script makes a plot of the CDF and PDF of t_coll for the ext1 model
# that is, the model nchain=9, resonance 5:4 (p=4)

# keep using the chainmodels module as a directory of metadata
import chainmodels
import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt
import scipy.stats

tkepler = chainmodels.tkepler

# This shoud be moved into a separate module, as several scripts can use it.
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

  lifetable = np.zeros(len(chainmodels.runs), 
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
filename = 'chain_data_library2.h5'
models, lifetable = build_lifetable(chainmodels.runs_ext1, filename)
filename = 'chain_data_library.h5'
models_old, lifetable_old = build_lifetable(chainmodels.runs, filename)

# merge the lifetables, but not the model list, so we will below only work with the one model 
#  from ext1
lifetable = np.vstack((lifetable, lifetable_old))


# this is a general form, here only used with one model
tlifes = {}
nfinished = {}
nlong = {}
nunfinished = {}
for model in models:
  selection = lifetable[np.nonzero(
                 np.logical_and(
                 np.logical_and(lifetable[:]['nchain']==model[0], lifetable[:]['p']==model[1]),
                 lifetable[:]['tlife']>0.0
                 )
                 )]
  tlifes[model]  = selection[:]['tlife']
  nfinished[model] = len(selection)
  nlong[model]= np.count_nonzero(
                 np.logical_and(
                 np.logical_and(lifetable[:]['nchain']==model[0], lifetable[:]['p']==model[1]),
                 lifetable[:]['lasttime']>0.9999e8) )
  nunfinished[model] = np.count_nonzero(
                 np.logical_and(
                 np.logical_and(lifetable[:]['nchain'] == model[0], lifetable[:]['p'] == model[1]),
                 np.logical_and(lifetable[:]['lasttime'] < 0.9999e8, lifetable[:]['tlife'] < 0.0)
                 ) )



def mensqerrkde(data,bwidth):
  '''
  Defines the cross-validation penalty function to minimize for the KDE bandwidth
  '''
  kde = scipy.stats.gaussian_kde(data, bw_method = bwidth) 
  # compute integral( kde^2 ) - 2/N sum_{j=1}^{N} (kdem1(j))(data[j]) 
  # limits of np.inf don't work very well, misses the action
  ikde2, quaderr = scipy.integrate.quad((lambda x: kde(x)*kde(x)), 
                     0.5*data.min(), 2*data.max(), points=np.median(data))
  # this loop does leave-one-out KDE and gets the value at the left-out location
  # np.roll shifts the array, and we always leave out the first data point
  kdesum = 0.0
  for j in range(0, len(data)):
    kdem1 = scipy.stats.gaussian_kde(data[1:], bw_method = bwidth) 
    kdesum += kdem1(data[0])[0]
    data = np.roll(data, shift=1)
  estimate_msqe = ikde2 - 2.0/(len(data)) * kdesum
  return estimate_msqe


# Now plot the CDF, and do the KDE based PDF estimation
for model in models:
  cdf_samples = np.linspace(0.0, 1.0, num = len(tlifes[model]), endpoint = False )
  cdf_samples += cdf_samples[1]
  
  plt.plot(np.sort(tlifes[model]), cdf_samples)
  plt.axvline(x = np.median(tlifes[model]), dashes=(5,5), color='grey')
  plt.title('Empirical CDF for mode '+ str(model))
  plt.ylabel('Emperical CDF')
  plt.xlabel(r'$t_{\rm coll}$')
  plt.ylim((0.0, 1.0))

  x = np.linspace(0.01,0.2)
  y = np.zeros(x.shape)
  for i,p in enumerate(x):
    y[i] = mensqerrkde(tlifes[model], x[i])

  plt.figure()
  plt.plot(x,y)
  plt.title('Cross-validation test for bandwidth factor')
  plt.xlabel('bandwidth factor')
  plt.ylabel('penalty')
  
  optimal_bandwidth = x[np.argmin(y)]
  print('Found optimal bandwidth ', optimal_bandwidth)
  plt.axvline(x=optimal_bandwidth, color='grey')

  plt.figure()
  trange = np.linspace(0.0*tlifes[model].min(), 1.1*tlifes[model].max(), num=500)
  kde = scipy.stats.gaussian_kde(tlifes[model], bw_method = optimal_bandwidth) 
  plt.plot(trange, kde(trange))
  plt.plot(tlifes[model], 0.0*tlifes[model], 'o')
  plt.axvline(x = np.median(tlifes[model]), dashes=(5,5), color='grey')
  plt.ylim(bottom=0)
  plt.xlim(left=0)
  plt.title('KDE estimate of lifetime PDF')
  plt.xlabel(r'$t_{\rm coll}$')
  plt.ylabel('Estimated PDF')

plt.show()
