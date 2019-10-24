# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign001 check the history of a given model
# 
import sys
import time
import json
import numpy as np
import rebound
import chaincalc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['axes.formatter.limits'] = (-2, 2)

from campaign001 import *

#Now have name col from campaign001
col = Campaign001()
print('Campaign has {} models'.format(col.get_size()))

nchain = int(sys.argv[1])
p_res = int(sys.argv[2])
seq = int(sys.argv[3])

print('looking for nchain={} p_res={} seq={}'.format(nchain, p_res, seq))

hh = col.indexes['p_res'][p_res]['nchain'][nchain]['seq'][seq]
print('Model hash ',hh)

m = col.get_model_hash(hh)

print('restoring Rebound from ',m.simarchive_filename)

sa = rebound.SimulationArchive(m.simarchive_filename)

print('Number of snapshots {:d}'.format(len(sa)))
print('Time limits {:e} {:e}'.format(sa.tmin, sa.tmax))
#print(dir(sa), sa.t[0:sa.nblobs])

a = np.zeros([nchain,len(sa)])
P = np.zeros([nchain,len(sa)])
l = np.zeros([nchain,len(sa)])
Omega = np.zeros([nchain,len(sa)])
omega = np.zeros([nchain,len(sa)])

for ip in range(0,nchain):
    for i, sim in enumerate(sa):
        a[ip,i] = sim.particles[ip+1].a     
        P[ip,i] = sim.particles[ip+1].P     
        l[ip,i] = sim.particles[ip+1].l     
        omega[ip,i] = sim.particles[ip+1].omega
        Omega[ip,i] = sim.particles[ip+1].Omega 

for ip in range(0,nchain-1):
   plt.plot(np.array(sa.t[0:sa.nblobs])/keplertime, a[ip,:], '-', label='{}'.format(ip+1))
plt.xlim(left=0)
plt.title('semimajor axis')
plt.legend()
    
plt.figure()
for ip in range(0,nchain-1):
   plt.plot(np.array(sa.t[0:sa.nblobs])/keplertime, P[ip+1,:]/P[ip+0,:],'-', label='{}:{}'.format(ip+1+1, ip+1))

plt.axhline(y=(p_res-1+1.0)/(p_res-1),color='grey',dashes=(4,4),linewidth=1)
plt.axhline(y=(p_res+1+1.0)/(p_res+1),color='grey',dashes=(4,4),linewidth=1)
plt.axhline(y=(p_res+1.0)/p_res,color='grey',dashes=(4,4),linewidth=1)
plt.xlim(left=0)
plt.title('Period ratios')
plt.legend()

def wrap2pi(x):
    return np.mod(x, 2.0*np.pi) - np.pi

def compresangles(l,Omega,omega,p,inner,outer):
    #varphi = Omega + omega
    outp = p+1
    theta1 = wrap2pi( outp*l[outer,:] -p*l[inner,:] 
                      - (Omega[outer,:] +omega[outer,:]) )
    theta2 = wrap2pi( outp*l[outer,:] -p*l[inner,:] 
                    - (Omega[inner,:] +omega[inner,:]) )
    return [theta1, theta2]


print('tdep ', m.params['tdep']/keplertime, 'deltatdep', m.params['deltatdep']/keplertime)
for ip in range(0,nchain-1):
    angles = compresangles(l, Omega, omega, p_res, ip, ip+1)

    plt.figure()
    plt.title('resonant angles {} and {}'.format(ip+1+1, ip+1))
    plt.scatter(np.array(sa.t[0:sa.nblobs])/keplertime, angles[0], s=1)
    plt.scatter(np.array(sa.t[0:sa.nblobs])/keplertime, angles[1], s=1)
    plt.ylim([-np.pi, np.pi])
    plt.xlim(left=0)
    plt.axvline(x=m.params['tdep']/keplertime, dashes=(4,4), color='grey')
    plt.axvline(x=(m.params['tdep']+m.params['deltatdep'])/keplertime, dashes=(4,4), color='grey')

plt.show()
