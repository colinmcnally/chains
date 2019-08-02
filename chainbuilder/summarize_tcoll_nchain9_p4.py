# This script makes a plot of the CDF and PDF of t_coll for the ext1 model
# that is, the model nchain=9, resonance 5:4 (p=4)
#
# We show that with a cut on the right (long time) tail fo the data, the remaining core is
# pretty indistinguishable from lognormal.
# 2008EP&S...60..681E have claimed that for random non-resoant systems with constant semimajor axis distributions
# the orbit crossing time distribution is lognormal, with a "dispersion" of 0.2. I think they mean a "scale" parameter
# relative to a lognormal with standard deviation unity. In that case, this data set has a lognormal core with 
# scale parameter 0.156.
#
#@ARTICLE{2008EP&S...60..681E,
#       author = {{Emori}, H. and {Nakazawa}, K. and {Iwasaki}, K.},
#        title = "{Probability distribution of orbital crossing times in a protoplanetary system}",
#      journal = {Earth, Planets, and Space},
#         year = "2008",
#        month = "Jun",
#       volume = {60},
#        pages = {681-691},
#       adsurl = {https://ui.adsabs.harvard.edu/abs/2008EP&S...60..681E},
#      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
#}

# keep using the chainmodels module as a directory of metadata
import chainmodels
import numpy as np
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


# Now plot the CDF, and do the KDE based PDF estimation
for model in models:
  cdf_samples = np.linspace(0.0, 1.0, num = len(tlifes[model]), endpoint = False )
  cdf_samples += cdf_samples[1]

  # transform data to log space 
  tlifes[model] = np.log(tlifes[model])
  
  plt.step(np.sort(tlifes[model]), cdf_samples, where='post')
  plt.axvline(x = np.median(tlifes[model]), dashes=(5,5), color='grey')
  plt.title('Empirical CDF for mode '+ str(model))
  plt.ylabel('Emperical CDF')
  plt.xlabel(r'$\ln(t_{\rm cross})$')
  plt.ylim((0.0, 1.0))

  grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                    {'bandwidth': np.logspace(-2, 0, 300)},
                    cv=LeaveOneOut(), iid=False) # cv-fold cross-validation
  grid.fit(tlifes[model].reshape(-1,1))
  print('GridSearchCV best params found ',grid.best_params_)

  trange = np.linspace(0.1*tlifes[model].min(), tlifes[model].max(), num=5000)

  kdesk = KernelDensity(kernel='gaussian', 
                        bandwidth=grid.best_params_['bandwidth'] ).fit(tlifes[model].reshape(-1,1) )
  integrated_kde_probability = scipy.integrate.trapz( 
           np.exp(kdesk.score_samples( trange.reshape(-1,1) ))*(1.0/np.exp(trange)),
           x= np.exp(trange)) 
  print('sklearn+gaussian total probability check', integrated_kde_probability )

  quartiles = np.percentile(tlifes[model], (25.0, 50.0, 75.0))
  cutindex = -10 # just lop off the last 10 visually.
  print('cut index ', cutindex,' of len ',len(tlifes[model]))

  # fit an elliptic envelope assuming the core is Gaussian
  # This is pointless, as it just identifies the same tail, and the level is set by the contammination param.
  cov = EllipticEnvelope(contamination=0.15, random_state=12).fit(tlifes[model].reshape(-1,1))
  covcut = np.argmin( cov.predict(np.sort(tlifes[model]).reshape(-1,1)) > 0.0)
  print('EllipticEnvelope cut',covcut)

  # fit a Minimum Covariance Determinant (MCD) robust estimator to data
  robust_cov = MinCovDet(random_state=12).fit(tlifes[model].reshape(-1,1))
  print('MinCovDet estimate of std dev',np.sqrt(robust_cov.covariance_[0][0]))
  
  normfit = scipy.stats.norm.fit(np.sort(tlifes[model])[:cutindex])
  print('MLE Gaussian fit parameters',normfit)
  norm = scipy.stats.norm(loc=normfit[0], scale= normfit[1])

  plt.figure()
  plt.plot(np.exp(trange), np.exp(kdesk.score_samples( trange.reshape(-1,1) ))*(1.0/np.exp(trange)), 
           label='sklearn+gaussian')
  plt.plot(np.exp(trange), norm.pdf( trange )*(1.0/np.exp(trange)), 
           label='MLE lognorm -tail')
  plt.plot(np.exp(tlifes[model]), 0.0*tlifes[model], '+')
  plt.axvline(x = np.exp(np.median(tlifes[model])), dashes=(5,5), color='grey')
  plt.ylim(bottom=-1e-8)
  plt.xlim(left=0)
  plt.title('KDE estimate of lifetime PDF for n=9  5:4 MMR')
  plt.xlabel(r'$t_{\rm cross}$')
  plt.ylabel('Estimated PDF')
  plt.legend()


  plt.figure()
  plt.axvline(x = np.median(np.exp(tlifes[model])), dashes=(5,5), color='grey')
  plt.step(np.sort(np.exp(tlifes[model])), cdf_samples, where='post', label='observations')
  plt.plot(np.exp(trange), 
           scipy.integrate.cumtrapz(np.exp(kdesk.score_samples( trange.reshape(-1,1) ))*(1.0/np.exp(trange)), 
           x= np.exp(trange), initial=0.0) , label = 'KDE in log space')
  plt.ylim((0.0, 1.0))
  plt.xlim(left=0.0)
  plt.title('CDF for mode '+ str(model))
  plt.ylabel('CDF')
  plt.xlabel(r'$t_{\rm cross}$')
  plt.legend()

  #Do some frequentist tests to see if the data is consistent with a couple distributions 
  # Forst, chekc the core fo the distribution is Gaussian
  astests = {}
  astests['norm'] = scipy.stats.anderson(
                      np.sort(tlifes[model]-np.median(tlifes[model]))[:cutindex], dist = 'norm')
  print(astests['norm'])
  print(' Anderson-Darling test for a lognormal distribution, discarding right tail ')
  print(' Test statistic {:3.2e}'.format(astests['norm'][0]))
  critindex = -1
  for i, crit in enumerate(astests['norm'][1]):
    if astests['norm'][0] > crit:
      critindex = i
  if critindex < 0:
    print(' Null hypothesis (samples are lognormal) accepted at {:3.1f}%'.format(astests['norm'][2][0]))
  else:
    print(' Null hypothesis (samples are lognormal) rejected at {:3.1f}%'.format(
                                                                         astests['norm'][2][critindex]))

  # test if full linear space data is Gumbel
  print(' ')
  #astests['gumbel'] = scipy.stats.anderson(np.exp(tlifes[model]), dist = 'gumbel')
  astests['gumbel'] = scipy.stats.anderson((tlifes[model]), dist = 'gumbel')
  print(astests['gumbel'])
  if astests['gumbel'][0] > astests['norm'][1][-1]:
    print(' Null hypothesis (samples are Gumbel distributed) rejected at {:3.1f}%'.format(
                                                                         astests['norm'][2][-1]))
  
  # make QQ plots of the full data and data with the tail cut off
  fig = plt.figure(figsize=(8,4))
  ax = [plt.subplot(121), plt.subplot(122)]
  scipy.stats.probplot(np.sort(tlifes[model]), dist='norm', plot=ax[0]) 
  ax[0].set_title('Normal Quantile-Quantile plot -- all data ')
  ax[0].set_ylabel(r'Observed $\log(t_{\rm cross})$')
  scipy.stats.probplot(np.sort(tlifes[model])[:cutindex], dist=norm, plot=ax[1]) 
  ax[1].set_title('Normal Quantile-Quantile plot discarding right tail')
  ax[1].set_ylabel(r'Observed $\log(t_{\rm cross})$')
  plt.tight_layout()

plt.show()
