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
tlifes = {}
for model in models:
  selection = lifetable[np.nonzero(
                 np.logical_and(
                 np.logical_and(lifetable[:]['nchain']==model[0], lifetable[:]['p']==model[1]),
                 lifetable[:]['tlife']>0.0
                 )
                 )]
  var[model]     = np.var(selection[:]['tlife'])
  tlifes[model]  = selection[:]['tlife']
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
plt.colorbar(label='mean(log( time to collison ))')
plt.title('Time to collision in resonant chains with q=$10^{-5}$')

plt.tight_layout()

fig = plt.figure(figsize=(8,7))
for i,model in enumerate(models):
  ax = fig.add_subplot(4,8,(x[i]-3)+8*(3-(y[i]-3))+1, xticks=[],yticks=[1e5, 1e6, 1e7, 1e8])
  ax.set_yscale('log')
  vals = np.sort(tlifes[model])
  if len(vals) > 0 :
    ax.semilogy(vals,'-o',markersize=2, color='C0')
  if len(vals) < 10 :
    ax.semilogy(np.arange(len(vals),10), [1e8] *(10-len(vals)),'^',markersize=4,color='C1')
  ax.set_ylim((5e4,2e8))
  if x[i] > 3:
    ax.set_yticklabels([])
  else:
    ax.set_ylabel('{:d}:{:d}  '.format(y[i]+1, y[i]) +r'$t_{\rm cross}$')
  if y[i] == 3:
    ax.set_xlabel(r'$n_{\rm chain} =$'+'{:d}'.format(x[i]))
fig.suptitle('Time to orbit crossing in resonant chains with q=$10^{-5}$')
fig.subplots_adjust(wspace=0.05, hspace = 0.05, top=0.93, bottom=0.05, right=0.95, left=0.1)
#plt.tight_layout()

plt.show()
