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
     'tau_e':-1e3,
     'tau_inc':-1e3,
     'a0':0.11,
     'redge':0.1,
     'deltaredge':0.1,
     'tdep':1e-10,
     'deltatdep':1e-10,
     'pmass':1e-5,
     'nchain':2,
     'q_res':1,
     'p_res':1,
     'seq':0,
     'collision':'line',
     'G':1.0,
     'starmass':1.0,
     'integrator':'WHFAST',
     'integrator_dt':1e-2*2*np.pi*0.1**1.5,
     'snap_wall_interval':15,
     'incscatter':np.pi/100.0,
     'aspread':0.05 }
p['physical_outputs'] = [p['tdep'], p['tdep']+p['deltatdep']]

def runtest(p):
    c = chaincalc.TauDampModel(p)

    c.init_rebound()

    print(c.hash)
    print(c.sim.status())

    c.set_wall_start(wall_start)
    try:
        for i in range(0,50000):
            c.sim.integrate((i+1)*100.0)
            c.manual_snapshot()
    except rebound.Collision:
        c.mark_collision()
    else:
        c.manual_snapshot()
    print('')

print('Running test without collision')
p['tdep'] = 1e8
p['deltatdep'] = 1e3
runtest(p)



