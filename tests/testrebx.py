# A prototype test of using reboundx, to help work out what the chain model class will do
import numpy as np
import sys
sys.path.append('/Users/colinm/vc/rebound')
sys.path.append('/Users/colinm/vc/reboundx')
import rebound
import reboundx
import time

wall_start = time.time()

keplertime = 2.0 * np.pi * 0.1**1.5

sim = rebound.Simulation() 
sim.add(m=1., hash="sun")

def add_planet(sim, m, a, true_anomaly, hash):
    # vary f as the true anomaly 
    sim.add(rebound.Particle(simulation=sim, primary=sim.particles[0], m=m, 
          a=a, e=0.0, inc = 0.0, Omega=0.0, omega=0.0, f = true_anomaly ,
          r = a*np.sqrt(m/3.0), hash=hash))

add_planet(sim, 1e-5, 1.0, 0.0, "earth")
add_planet(sim, 1e-5, 1.05, np.pi, "e2")

#sim.collision = "line"
sim.collision = "direct"

# To be able to reference something from the callback, pass the method on some object, then we get self
# available inside the method!
# Do all the I/O from inside this object, so we can just call 
class simrecorder:

    def __init__(self, wall_start, heart_print_interval=5.0):#, snap_wall_interval=10.0):
        self.heart_print_interval = heart_print_interval
        self.lastheart = wall_start
        #self.snap_wall_interval = snap_wall_interval
        #self.lastsnap = wall_start
        self.wall_start = wall_start

    def heartbeat(self,sim):
        # sim is a pointer, so use contents 
        wall_now = time.time()
        if (wall_now - self.lastheart > self.heart_print_interval): 
            self.lastheart = wall_now
            print("time {} walltime {}".format(sim.contents.t, wall_now -self.wall_start), end='\r')
        #if (wall_now - self.lastsnap > self.snap_wall_interval):
        #    self.lastsnap = wall_now
        #    print("\nsnapshot at {}".format(wall_now -self.wall_start))

rebx = reboundx.Extras(sim)
mof = rebx.load_force("modify_orbits_resonance_relax") 
mof.params['tau_a'] = 0.035
mof.params['res_aspectratio0'] = 0.035
mof.params['res_sigma0']       = 3.8e-4 *0.5*1.4
mof.params['res_redge']        = 0.1
mof.params['res_deltaredge']   = 0.001
mof.params['res_alpha']        = 1.0
mof.params['res_flaringindex'] = 2.0/7.0
mof.params['res_ffudge']       = 1.0/100.0
mof.params['res_tdep']         = 4e5 * keplertime
mof.params['res_deltatdep']    = 1e3 * keplertime
rebx.add_force(mof)

snap_wall_interval = 10.0
simarchive_filename = "simulationarchive.rbsa"

ht = simrecorder(wall_start)
sim.heartbeat = ht.heartbeat
sim.automateSimulationArchive(simarchive_filename, walltime=snap_wall_interval, deletefile=True)


try:
  sim.integrate(100000.0)
except rebound.Collision:
    collided = []
    for p in sim.particles:
        if p.lastcollision == sim.t:
            collided.append(p.index)
    # Custom resolution
    print("") # newline
    print("Particles {} collided sim.t {}".format(collided, sim.t))
    sim.simulationarchive_snapshot(simarchive_filename)
else:
    sim.simulationarchive_snapshot(simarchive_filename)


print("") # newline
print('status',sim.status())

sa = rebound.SimulationArchive(simarchive_filename)
print("Number of snapshots {}".format(len(sa)))
print("Time of first {:e} last {:e}".format(sa.tmin, sa.tmax))
