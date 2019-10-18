# Colin McNally 2019 <colin@colinmcnally.ca>
# 
# A class that defines a resonant chain model
#
# Needs to define:
#  - parameters for the run
#  - a unique hash used to lookup the model, identify files, etc. Need to be able to calculate this without setting up the rest.
#  - the rebound  simulation object with the parameters and function pointers set
#  - the reboundx object with the forcing parameters set
#  - the heartbeat callback
#  - the simulation integrator method, including checkpoint, releading, and status I/O to disc
#  - a method to check the simulation status from disc
#
#  A given set of modelsi (i.e. their hashes) will constitute a simulation campaign, and the driver will be a
#  class which figure out how to get all the models to be completed
import sys
import os.path
import hashlib
import time
import random
import json
import numpy as np
import rebound
import reboundx

class Model:
    def __init__(self, params, snap_wall_interval=60*60):
        """ Params is a dict of parameters, easier to reuse when initializing in a loop
            aspectratio0
            sigma0
            redge
            deltaredge
            alpha
            flaringindex
            ffudge
            tdep
            deltatdep
            pmass
            nchain
            p_res
            a0
            seq

            some other values set as object attributes, but not intended to be modified
            G
            starmass
            integrator
            integrator_dt
     """
        paramlist = ['aspectratio0',
                     'sigma0',
                     'redge',
                     'deltaredge',
                     'alpha',
                     'flaringindex',
                     'ffudge',
                     'tdep',
                     'deltatdep',
                     'pmass',
                     'nchain',
                     'p_res',
                     'q_res',
                     'a0',
                     'seq']
                     
        for parstr in paramlist:
            try:
                parval = params[parstr]
            except KeyError as keyerr:
                print('chainmodel __init__ did not find param {}'.format(keyerr.args[0]))
                raise
            else:
                setattr(self, parstr, parval)
        #probably keep the above setattr, and remove this, or do the opposite
        self.params = params
        self.G = 1.0
        self.starmass = 1.0
        self.integrator = "whfast" 
        # 1e-2 times the inner edge orbital time
        self.integrator_dt = 1e-2*2.0*np.pi* self.redge**1.5
        self.set_model_hash()
        self.simarchive_filename = os.path.join('output','sim'+self.hash+'.rbsa')
        self.status_filename = os.path.join('output','sim'+self.hash+'_status.json')

        self.heart_print_interval = 5.0
        self.snap_wall_interval = snap_wall_interval

    def set_wall_start(self, wall_start):
        self.wall_start = wall_start
        self.lastheart = wall_start

    def set_model_hash(self):
        """Calculate a hash specifying the model. Need to use all the parameters of the model here."""
        params = self.params
        repstr = ''
        for pk in sorted(params.keys()):
            repstr += '{:e}'.format(params[pk])
        repstr += '{:e}'.format(self.G)
        repstr += '{:e}'.format(self.starmass)
        repstr += '{}'.format(self.integrator)
        repstr += '{:e}'.format(self.integrator_dt)
        #calculate the MD5 hash of this string
        self.hash = hashlib.md5(repstr.encode(encoding='utf-8')).hexdigest()


    def lock(self):
        """ Set lock in status file """
        self.status['lock'] = True
        self.overwrite_status()

    def unlock(self):
        """ Unset lock in status file """
        self.status['lock'] = False
        self.overwrite_status()

    def overwrite_status(self):
        """Overwrite the status file, it's a JSON format dict"""
        self.status['sim.t'] = self.sim.t
        # formally, the params are redundant, as you can;t calculate the filename hash
        #  without knowing them!
        self.status['params'] = self.params
        # Could add version, time, and host info here
        with open(self.status_filename,'w') as json_file:
            json.dump(self.status, json_file)
        
    def get_status(self):
        """Try to read the status JSON file, fail and pass through FileNotFound"""
        try:
            with open(self.status_filename,'r') as json_file:
                self.status = json.load(json_file)
                print('loaded status')
        except FileNotFoundError:
            print('No old status file was found at: {}'.format(self.status_filename))
            raise
  
    def init_rebound(self):
        """Initialize the rebound simulation and rebx object for this model """
        #in future, this will take an argument specifying to try to reload from disc

        try:
            self.get_status()
            #Here, it is possible that self.sim.t is > status['sim.t] is the last achive was done automatically,
            # but that will nto cuase a problem in current usage. It might be better to move all checkpointing
            # to the driver integrate loop where we can control everything, as there's no callback in REBOUND.
            self.sim = rebound.Simulation(self.simarchive_filename) 
            self.sim.automateSimulationArchive(self.simarchive_filename, walltime=self.snap_wall_interval, deletefile=False)
        except FileNotFoundError:
            self.sim = rebound.Simulation()
            self.sim.integrator = self.integrator
            self.sim.collision = "line"
            self.sim.automateSimulationArchive(self.simarchive_filename, walltime=self.snap_wall_interval, deletefile=True)
            # status should be running or collided
            self.status = {'status':'running'}
            self.unlock() # calls overwrite_status, and allows locking by driver

            # here we do the initialization for a chain model, maybe make this a method so we can overload only it
            self.sim.add(m=self.starmass, hash='star')

            #init random number generator to get reproducible random phases
            rng = random.Random(self.hash)
            a = self.a0
            for ip in range(1, self.nchain+2):
                true_anomaly = rng.uniform(0, 2.0*np.pi)
                print('adding a={:e}'.format(a))
                pt = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],
                                      m=self.pmass, a=a, f=true_anomaly)
                pt.r = a*np.sqrt(pt.m/3.0)
                self.sim.add(pt)
                # the 0.05 is a spread. but maybe this could be a param
                a *= (self.p_res / (self.p_res + self.q_res + 0.05))**(-2.0/3.0)
         
        #always reset function pointers
        self.sim.heartbeat = self.heartbeat 

        self.rebx = reboundx.Extras(self.sim)
        mof = self.rebx.load_force("modify_orbits_resonance_relax")
        mof.params['res_aspectratio0'] = self.aspectratio0
        mof.params['res_sigma0'] = self.sigma0
        mof.params['res_redge'] = self.redge
        mof.params['res_deltaredge'] = self.deltaredge 
        mof.params['res_alpha'] = self.alpha
        mof.params['res_flaringindex'] = self.flaringindex
        mof.params['res_ffudge'] = self.ffudge
        mof.params['res_tdep'] = self.tdep
        mof.params['res_deltatdep'] = self.deltatdep
        self.rebx.add_force(mof)
        self.mof = mof

       

    def heartbeat(self, sim):
       """Heartbeat callback from rebound integrate
          sim is a pointer (see ctypes) not the actual object """
       wall_now = time.time()
       if (wall_now - self.lastheart > self.heart_print_interval):
            self.lastheart = wall_now
            print("time {} walltime {}".format(sim.contents.t, wall_now -self.wall_start), end='\r')


    def mark_collision(self):
        """Called when the driver gets a collision exception
           Sets the model status, and records the collided particles 
           Forces overwrite of the status file """
        collided = []
        for p in self.sim.particles:
            if p.lastcollision == self.sim.t:
                collided.append(p.index)
        print('') # newline
        print('Particles {} collided sim.t {}'.format(collided, self.sim.t))
        self.status['status'] = 'stop collision'
        self.status['collided'] = collided
        self.overwrite_status()
 

    def manual_snapshot(self):
        """Call to save a manual snapshot, for example when driver hits walltime limit"""
        self.sim.simulationarchive_snapshot(self.simarchive_filename)
        self.overwrite_status()
