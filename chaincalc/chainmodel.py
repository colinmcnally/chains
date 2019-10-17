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
import hashlib
import numpy as np

class model:
    def __init__(self, params):
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
                     'p_res']
                     
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

    def set_model_hash(self):
        """Calculate a hash specifying the model. Need to use all the parameters of the model here
        """
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
