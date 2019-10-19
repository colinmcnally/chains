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

#this is sort of a dummy driver

wall_start = time.time()

p = {'aspectratio0':0.035,
     'sigma0':3.8e-4,
     'redge':0.1,
     'a0':0.11,
     'deltaredge':0.001,
     'alpha':1.0,
     'flaringindex':2.0/7.0,
     'ffudge':1.0/100.0,
     'tdep':1e1,
     'deltatdep':1e1,
     'pmass':1e-5,
     'nchain':4,
     'q_res':3,
     'p_res':40,
     'seq':0 }

def runtest(p):
    c = chaincalc.Model(p)

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



