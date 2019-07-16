# calls custom Rebound setup, as a library, which outputs to a HDF5 file

import chainmodels
import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt


lifetimes = {}
for key in chainmodels.filenames:
  filename = chainmodels.filenames[key]
  #print(filename)
  orb = h5py.File(filename, 'r')
  tcoll = orb.attrs['tcollstop'][0]
  tstart = orb.attrs['tdep'][0]+orb.attrs['deltatdep'][0]
  if tcoll > 0.0:
    tlife = tcoll - tstart
  else:
    tlife = 1e8
  #print('tlife {:e}'.format(tlife))
  lifetimes[key] = tlife
  orb.close()

lifetable = np.zeros(len(chainmodels.runs), 
                     dtype=[('nchain',np.int16), ('p',np.int16), ('realization',np.int16), ('tlife',np.float64)])
for i,key in enumerate(chainmodels.runs):
  nchain, p, realization = chainmodels.runs[key]
  lifetable[i]['nchain'] = nchain
  lifetable[i]['p'] = p
  lifetable[i]['realization'] = realization
  lifetable[i]['tlife'] = lifetimes[key]


models = set()
for i,key in enumerate(chainmodels.runs):
  nchain, p, realization = chainmodels.runs[key]
  models.add((nchain,p))

var = {}
meanlog = {}
for model in models:
  selection = lifetable[np.nonzero(np.logical_and(lifetable[:]['nchain']==model[0], lifetable[:]['p']==model[1]))]
  var[model]     = np.var(selection[:]['tlife'])
  meanlog[model] = np.log10(selection[:]['tlife']).mean()
  

x = np.zeros(len(models), dtype=np.int16)
y = np.zeros(len(models), dtype=np.int16)
c = np.zeros(len(models))
s = np.zeros(len(models))
for i,model in enumerate(models):
  nchain, p = model 
  x[i] = nchain
  y[i] = p
  c[i] = meanlog[model]
  s[i] = 50.0*np.log10(var[model]+1) +20.0

plt.scatter( x, y, c=c, s=s)
plt.yticks([3,4,5,6], ['4:3','5:4','6:5','7:6'])
plt.xlabel('Number of planets in chain')
plt.ylabel('Resonance')
plt.colorbar(label='mean(log( lifetime))')
plt.title('Time to collision in resonant chains with q=$10^{-5}$')

plt.tight_layout()
plt.show()
