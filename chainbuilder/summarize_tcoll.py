# make tcoll plot from trajectory library file
# keep using the chainmodels module as a directory of metadata
import chainmodels
import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt



tkepler = 2*np.pi*0.1**1.5

library_file = h5py.File('chain_data_library.h5', 'r')

lifetimes = {}
lasttimes = {}
for key in chainmodels.filenames:
  filename = chainmodels.filenames[key]
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
                     dtype=[('nchain',np.int16), ('p',np.int16), ('realization',np.int16), ('tlife',np.float64), ('lasttime',np.float64)])
for i,key in enumerate(chainmodels.runs):
  nchain, p, realization = chainmodels.runs[key]
  lifetable[i]['nchain'] = nchain
  lifetable[i]['p'] = p
  lifetable[i]['realization'] = realization
  lifetable[i]['tlife'] = lifetimes[key]
  lifetable[i]['lasttime'] = lasttimes[key]


models = set()
for i,key in enumerate(chainmodels.runs):
  nchain, p, realization = chainmodels.runs[key]
  models.add((nchain,p))


# this effectively does a set of table queries. Maybe there's a better waty of doing a SQL "select where" on HDF5
var = {}
meanlog = {}
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
  var[model]     = np.var(selection[:]['tlife'])
  meanlog[model] = np.log10(selection[:]['tlife']).mean()
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

x = np.zeros(len(models), dtype=np.int16)
y = np.zeros(len(models), dtype=np.int16)
c = np.zeros(len(models))
s = np.zeros(len(models))
for i,model in enumerate(models):
  nchain, p = model 
  x[i] = nchain
  y[i] = p
  c[i] = meanlog[model]
  s[i] = 50.0*np.log10(var[model]+1)

plt.figure(figsize=(10,4))
plt.scatter( x, y, c='black', s=5)
plt.scatter( x, y, c=c, s=s)
for i,model in enumerate(models):
    plt.annotate('L{:d}, U{:d}, F{:d}'.format(nlong[model], nunfinished[model],
                  nfinished[model]), (x[i],y[i]-0.35), fontsize=8)
plt.yticks([3,4,5,6], ['4:3','5:4','6:5','7:6'])
plt.xlabel('Number of planets in chain')
plt.ylabel('Resonance')
plt.ylim((2.5,6.3))
plt.xlim((2.5,10.95))
plt.colorbar(label='mean(log( lifetime))')
plt.title('Time to collision in resonant chains with q=$10^{-5}$')

plt.tight_layout()


plt.show()
