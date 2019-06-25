import numpy as np
import numpy.ma as ma
import sys
import os
import os.path
import h5py

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
mpl.rc('text', usetex=True)
plt.style.use('seaborn-colorblind')


#orbtxt = np.genfromtxt('orbits.txt',names=['t','a','e','inc','Omega','omega','l','P','f'])
#note:
# l = longitude
#varphi = Omega + omega
orb = h5py.File('orbits.h5','r')

nplanets = orb.attrs['nchain'][0]
print(nplanets)

plt.figure()
for i in range(0, nplanets):
  #plt.plot(planets[0]['t'], planets[i]['a'])
  plt.plot(orb['t'], orb['a'][i,:])
plt.title('Semi-major axis')


plt.figure()
#  plt.plot(p['t'],p['P'])
for i in range(0, nplanets-1):
  plt.plot(orb['t'], orb['P'][i+1,:]/orb['P'][i,:])
for p in range(2,8):
  y = (p + 1.0) / p
  plt.axhline(y=y, color='black', linewidth=1)
  plt.annotate('{:d}:{:d}'.format(p+1,p), (0.0,y))
plt.title('Period ratios')

plt.figure()
for i in range(0,nplanets):
  plt.plot(orb['t'],orb['e'][i,:])
plt.figure()
for p in range(0,nplanets):
  plt.plot(orb['t'],orb['inc'][i,:])

def wrap2pi(x):
  return np.mod(x, 2.0*np.pi) - np.pi

def compresangles(orbits,p):
  #varphi = Omega + omega
  outp = p+1
  theta1 = wrap2pi( outp*orb['l'][1,:] -p*orbits['l'][0,:] - (orbits['Omega'][1,:] +orbits['omega'][1,:]) )
  theta2 = wrap2pi( outp*orb['l'][1,:] -p*orbits['l'][0,:] - (orbits['Omega'][0,:] +orbits['omega'][0,:]) )
  return [theta1, theta2]

angles = compresangles(orb, orb.attrs['p'][0])

plt.figure()
plt.scatter(orb['t'], angles[0], s=1)
plt.scatter(orb['t'], angles[1], s=1)
plt.axvline(x = orb.attrs['tdep'][0], color='grey')
plt.axvline(x = orb.attrs['tdep'][0]+orb.attrs['deltatdep'][0], color='grey')
plt.ylim([-np.pi, np.pi])
plt.xlim(left=0)
plt.show()
