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
#  A given set of models (i.e. their hashes) will constitute a simulation campaign, and the driver will be a
#  class which figure out how to get all the models to be completed
import sys
import os
import os.path
import hashlib
import time
import random
import copy
import json
import numpy as np
import rebound
import reboundx

class Model:
    def validate_paramlist(self, params):
        """Check the param list, designed to be overloaded."""
        paramlist =  ['aspectratio0',
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
                      'seq',
                      'G',
                      'collision',
                      'starmass',
                      'integrator',
                      'integrator_dt',
                      'snap_wall_interval',
                      'incscatter',
                      'aspread',
                      'physical_outputs']
        for parstr in paramlist:
            try:
                parval = params[parstr]
            except KeyError as keyerr:
                print('chainmodel did not find param {}'.format(keyerr.args[0]))
                raise
 

    def __init__(self, params):
        """Params is a dict of parameters, easier to reuse when initializing in a loop.
           Init broken up into methods to make overloading easier."""
        self.validate_paramlist(params)
        self.other_init(params)


    def other_init(self, params):
        """Initialization to be done after the parameters are checked."""
        #probably keep the above setattr, and remove this, or do the opposite
        self.params = copy.deepcopy(params) # need to do a deep copy as the dict is actually object refs
        self.set_model_hash()
        self.simarchive_filename = os.path.join('output','sim'+self.hash+'.rbsa')
        self.status_filename = os.path.join('output','sim'+self.hash+'_status.json')
        self.heart_print_interval = 5.0


    def set_wall_start(self, wall_start):
        """Store the wallclock start time of the driver process, not the Rebound simulation """
        self.wall_start = wall_start
        self.lastheart = wall_start


    def set_model_hash(self):
        """Calculate a hash specifying the model. Need to use all the parameters of the model here."""
        params = self.params
        repstr = ''
        # got to be careful to catch all cases here
        for pk in sorted(params.keys()):
            val = params[pk]
            if isinstance(val, str):
                repstr += '{}'.format(val)
            elif isinstance(val, int):
                repstr += '{:d}'.format(val)
            elif isinstance(val, float):
                repstr += '{:e}'.format(val)
            else:
                repstr += '{}'.format(val)
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
            #only do this from the driver side with manual trigger, as we have to hit the walltime limit anyways
            #self.sim.automateSimulationArchive(self.simarchive_filename, walltime=self.params['snap_wall_interval'], deletefile=False)
        except FileNotFoundError:
            self.sim = rebound.Simulation()
            if self.params['integrator']=='WHFASTUNSAFE11':
                self.sim.integrator = 'whfast'
                self.sim.ri_whfast.safe_mode = 0
                self.sim.ri_whfast.corrector = 11
            else:    
                self.sim.integrator = self.params['integrator']
            self.sim.collision = self.params['collision']
            #only do this from the driver side with manual trigger, as we have to hit the walltime limit anyways
            #self.sim.automateSimulationArchive(self.simarchive_filename, walltime=self.params['snap_wall_interval'], deletefile=True)
            # if not calling auto simarchive, need to at least clear out any old file
            if os.path.isfile(self.simarchive_filename):
                os.remove(self.simarchive_filename)
            # status should be running or collided
            self.status = {'status':'running'}
            self.unlock() # calls overwrite_status, and allows locking by driver

            # here we do the initialization for a chain model, maybe make this a method so we can overload only it
            self.sim.add(m=self.params['starmass'], hash='star')

            #init random number generator to get reproducible random phases
            rng = random.Random(self.hash)
            a = self.params['a0']
            for ip in range(1, self.params['nchain']+1):
                true_anomaly = rng.uniform(0.0, 2.0*np.pi)
                inclination = rng.normalvariate(0.0, self.params['incscatter'])
                print('adding {} at a={:e}'.format(ip,a))
                pt = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],
                                      m=self.params['pmass'], a=a, inc=inclination, f=true_anomaly)
                pt.r = a*np.sqrt(pt.m/3.0)
                self.sim.add(pt)
                # the 0.05 is a spread. but maybe this could be a param
                a *= (self.params['p_res'] / (self.params['p_res'] + self.params['q_res'] + self.params['aspread']))**(-2.0/3.0)
            self.sim.move_to_com()
            self.manual_snapshot()

        #the reboundx forces are initialized following params, and in a method so they can be overidden
        self.add_force()
         
        #always reset function pointers
        self.sim.heartbeat = self.heartbeat 


    def add_force(self):
        """Add reboundx forces. This is a method so it can be overloaded easily"""
        print(' adding modify_orbits_resonance_relax')
        self.rebx = reboundx.Extras(self.sim)
        mof = self.rebx.load_force('modify_orbits_resonance_relax')
        mof.params['res_aspectratio0'] = self.params['aspectratio0']
        mof.params['res_sigma0'] = self.params['sigma0']
        mof.params['res_redge'] = self.params['redge']
        mof.params['res_deltaredge'] = self.params['deltaredge']
        mof.params['res_alpha'] = self.params['alpha']
        mof.params['res_flaringindex'] = self.params['flaringindex']
        mof.params['res_ffudge'] = self.params['ffudge']
        mof.params['res_tdep'] = self.params['tdep']
        mof.params['res_deltatdep'] = self.params['deltatdep']
        self.rebx.add_force(mof)
        self.mof = mof


    def heartbeat(self, sim):
        """Heartbeat callback from rebound integrate
           sim is a pointer (see ctypes) not the actual object """
        wall_now = time.time()
        if (wall_now - self.lastheart > self.heart_print_interval):
            self.lastheart = wall_now
            print("time {:e} walltime {:e}".format(sim.contents.t, wall_now -self.wall_start), end='\r')


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


class TauDampModel(Model):
    """A model using constanat timescale damping, with an inner edge to the disc """

    def validate_paramlist(self, params):
        """Checks required params are in the model parameter dict"""
        paramlist =  [
                      'tau_a',
                      'tau_e',
                      'tau_inc',
                      'redge',
                      'deltaredge',
                      'tdep',
                      'deltatdep',
                      'pmass',
                      'nchain',
                      'p_res',
                      'q_res',
                      'a0',
                      'seq',
                      'G',
                      'collision',
                      'starmass',
                      'integrator',
                      'integrator_dt',
                      'snap_wall_interval',
                      'incscatter',
                      'aspread',
                      'physical_outputs']
        for parstr in paramlist:
            try:
                parval = params[parstr]
            except KeyError as keyerr:
                print('chainmodel TauDampModel did not find param {}'.format(keyerr.args[0]))
                raise


    def add_force(self):
        """Add tau damping force"""
        print('adding modify_orbits_forces_edge')
        self.rebx = reboundx.Extras(self.sim)
        mof = self.rebx.load_force('modify_orbits_forces_edge')
        mof.params['res_redge'] = self.params['redge']
        mof.params['res_deltaredge'] = self.params['deltaredge']
        mof.params['res_tdep'] = self.params['tdep']
        mof.params['res_deltatdep'] = self.params['deltatdep']
        self.rebx.add_force(mof)
        self.mof = mof
        for pt in self.sim.particles[1:]:
            pt.params['tau_a'] = self.params['tau_a']
            pt.params['tau_e'] = self.params['tau_e']
            pt.params['tau_inc'] = self.params['tau_inc']
