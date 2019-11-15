# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Vizualize test output
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

status_filename = sys.argv[1]

with open(status_filename,'r') as json_file:
  status = json.load(json_file)

sa = rebound.SimulationArchive(status['simarchive_filename'])

oa = chaincalc.OrbitArray(sa)

keplertime = 2.0*np.pi*1.0**1.5


print('Number of snapshots {:d}'.format(len(sa)))
print('Time limits {:e} {:e}'.format(sa.tmin, sa.tmax))

nchain = len(sa[0].particles) -1
print('nchain',nchain)


for ip in range(0,nchain):
   plt.plot(oa.orbitdata['t'][ip,:]/keplertime, oa.orbitdata['a'][ip,:],'-')
plt.xlim(left=0)
plt.title('semimajor axis')
    
plt.figure()
for ip in range(0,nchain-1):
   plt.plot(oa.orbitdata['t'][ip,:]/keplertime, oa.orbitdata['P'][ip+1,:]/oa.orbitdata['P'][ip,:],'-')
plt.xlim(left=0)
plt.title('Period ratios')

def wrap2pi(x):
    return np.mod(x + np.pi, 2.0*np.pi) - np.pi

def compresangles_old(l,Omega,omega,p,inner,outer):
    #varphi = Omega + omega
    outp = p+1
    theta1 = wrap2pi( outp*l[outer,:] -p*l[inner,:] 
                      - (Omega[outer,:] +omega[outer,:]) )
    theta2 = wrap2pi( outp*l[outer,:] -p*l[inner,:] 
                    - (Omega[inner,:] +omega[inner,:]) )
    return [theta1, theta2]

def compresangles(inner, outer, p):
    #varphi = Omega + omega
    outp = p+1
    theta1 = wrap2pi( outp*outer['l'] -p*inner['l'] 
                      - (outer['Omega'] +outer['omega']) )
    theta2 = wrap2pi( outp*outer['l'] -p*inner['l'] 
                    - (inner['Omega'] +inner['omega']) )
    return [theta1, theta2]


# try detecting the resonance
# average only the last period 
tlook = 1e4*keplertime #probably 1e2-1e3 times the output cadence
lko = oa.tail_orbitarray(tlook)

# could be ratio of means, or mean of ratios...
pratios = lko['P'][1:,:].mean(axis=1) / lko['P'][:-1,:].mean(axis=1)
meanpratios = (lko['P'][1:,:] / lko['P'][:-1,:]).mean(axis=1)
varpratios = (lko['P'][1:,:] / lko['P'][:-1,:]).var(axis=1)
print(' mean',pratios,'meanpratios',meanpratios,'varpratios',varpratios)



finding = oa.is_this_a_constant_ratio_chain(oa.compute_tight_angles(lko))
if finding[0]:
    p_res = finding[1] 
    print("found constant period ratio chain {}".format(p_res))
    for ip in range(0,nchain-1):
        #angles = compresangles_old(lookorbit['l'], lookorbit['Omega'], lookorbit['omega'], , ip, ip+1)
        angles = oa.compresangles(lko[:][ip,:], lko[:][ip+1,:], p_res)

        plt.figure()
        plt.title('resonant angles {} and {}'.format(ip+1+1, ip+1))
        plt.scatter(lko['t'][ip,:]/keplertime, angles[0], s=1)
        plt.scatter(lko['t'][ip,:]/keplertime, angles[1], s=1)
        plt.ylim([-np.pi, np.pi])


plt.figure()
for ip in range(0,nchain-1):
  plt.scatter(oa.orbitdata['t'][ip,:]/keplertime, wrap2pi( oa.orbitdata['pomega'][ip+1]-oa.orbitdata['pomega'][ip]), s=1) 

plt.show()
