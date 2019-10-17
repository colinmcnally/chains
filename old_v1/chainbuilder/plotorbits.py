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
orb = h5py.File(sys.argv[1], 'r')

nplanets = orb.attrs['nchain'][0]
ie = orb.attrs['lastout'][0]
print(nplanets)

plt.figure()
for i in range(0, nplanets):
  #plt.plot(planets[0]['t'], planets[i]['a'])
  plt.plot(orb['t'][:ie], orb['a'][i,:ie])
  plt.axvline(x = orb.attrs['tdep'][0], color='grey')
  plt.axvline(x = orb.attrs['tdep'][0]+orb.attrs['deltatdep'][0], color='grey')
plt.title('Semi-major axis')


plt.figure()
#  plt.plot(p['t'],p['P'])
for i in range(0, nplanets-1):
  plt.plot(orb['t'][:ie], orb['P'][i+1,:ie]/orb['P'][i,:ie])
for p in range(2,8):
  y = (p + 1.0) / p
  plt.axhline(y=y, color='black', linewidth=1)
  plt.annotate('{:d}:{:d}'.format(p+1,p), (0.0,y))
plt.title('Period ratios')

plt.figure()
for i in range(0,nplanets):
  plt.plot(orb['t'][:ie],orb['e'][i,:ie])
plt.figure()
for p in range(0,nplanets):
  plt.plot(orb['t'][:ie],orb['inc'][i,:ie])

def wrap2pi(x):
  return np.mod(x, 2.0*np.pi) - np.pi

def compresangles(orbits,p,inner,outer):
  #varphi = Omega + omega
  outp = p+1
  theta1 = wrap2pi( outp*orb['l'][outer,:ie] -p*orbits['l'][inner,:ie] - (orbits['Omega'][outer,:ie] +orbits['omega'][outer,:ie]) )
  theta2 = wrap2pi( outp*orb['l'][outer,:ie] -p*orbits['l'][inner,:ie] - (orbits['Omega'][inner,:ie] +orbits['omega'][inner,:ie]) )
  return [theta1, theta2]



for inner in range(0,nplanets-1):
  angles = compresangles(orb, orb.attrs['p'][0], inner, inner+1)

  plt.figure()
  plt.scatter(orb['t'][:ie], angles[0], s=1)
  plt.scatter(orb['t'][:ie], angles[1], s=1)
  plt.axvline(x = orb.attrs['tdep'][0], color='grey')
  plt.axvline(x = orb.attrs['tdep'][0]+orb.attrs['deltatdep'][0], color='grey')
  plt.ylim([-np.pi, np.pi])
  plt.xlim(left=0)
  plt.title('{:d} and {:d} in {:d}:{:d}'.format(inner, inner+1, orb.attrs['p'][0], orb.attrs['p'][0]+1))

plt.show()
