# Colin McNally 2019 <colin@colinmcnally.ca>
# Script to test the Model class.
# Set to run in laptop, not container
# 
import time
import sys
sys.path.append('/Users/colinm/vc/rebound')
sys.path.append('../crbx')
sys.path.append('../')
import chaincalc
import rebound
import numpy as np

#this is sort of a dummy driver

wall_start = time.time()

p = {'tau_a':-1e5,
     'tau_e':-1e5/50.0,
     'tau_inc':-1e5/50.0,
     'redge':1.0,
     'set_a':1.0,
     'deltaredge':1.0,
     'tdep':1e5*2*np.pi,
     'deltatdep':1e4*2*np.pi,
     'pmass':1e-5,
     'nchain':6,
     'a':None,
     'seq':0,
     'collision':'direct',
     'G':1.0,
     'starmass':1.0,
     'integrator':'WHFAST',
     'integrator_dt':1e-2*2*np.pi*1.0**1.5,
     'snap_wall_interval':60*15,
     'incscatter':0.0,
     #'aspread':0.05,
     'physical_output_dt':100.0*2*np.pi }
p['physical_outputs'] = []

# set up chain in equal mass
a = [1.0]
for i in range(1,p['nchain']):
    a.append(a[i-1]*(6./5.)**(2./3.)*1.02)
p['a'] = a

def runtest(p):
    c = chaincalc.SetAModel(p, verbose=True)

    c.init_rebound()

    print(c.hash)
    print(c.sim.status())

    c.set_wall_start(wall_start)
    try:
        for i in np.arange(0.0,1.3e5*2*np.pi,p['physical_output_dt']):
            c.sim.integrate(i)
            c.manual_snapshot()
    except rebound.Collision:
        c.mark_collision()
    else:
        c.manual_snapshot()
    print('')

print('Running test without collision')
runtest(p)



