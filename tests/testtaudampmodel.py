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
     'deltaredge':1.0/100.0,
     'tdep':1e1,
     'deltatdep':1e1,
     'pmass':1e-5,
     'nchain':4,
     'q_res':3,
     'p_res':40,
     'seq':0,
     'collision':'line',
     'G':1.0,
     'starmass':1.0,
     'integrator':'WHFAST',
     'integrator_dt':1e-2*2*np.pi*0.1**1.5,
     'snap_wall_interval':60*60,
     'incscatter':np.pi/100.0,
     'aspread':0.05 }

def runtest(p):
    c = chaincalc.TauDampModel(p)

    c.init_rebound()

    print(c.hash)
    print(c.sim.status())

    c.set_wall_start(wall_start)
    try:
        c.sim.integrate(10000.0)
    except rebound.Collision:
        c.mark_collision()
    else:
        c.manual_snapshot()
    print('')

print('Running test with collision')
runtest(p)

print('Running test without collision')
p['p_res'] = 5
p['q_res'] = 1
p['tdep'] = 1e5
p['deltatdep'] = 1e3
runtest(p)



