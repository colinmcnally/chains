# Colin McNally 2019 <colin@colinmcnally.ca>
# Campaign008 definitions
# 
import numpy as np
import random
import sys
import chaincalc

keplertime = 2.0*np.pi*1.0**1.5
wall_check_interval = 15*60

class Campaign008(chaincalc.CampaignBase):
    def __init__(self):
        p = {'tau_a':-1e5,
             'tau_e':-1e5/50.0,
             'tau_inc':-1e5/50.0,
             'redge':1.0,
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
             'physical_output_dt':100.0*2*np.pi }
        p['physical_outputs'] = []

        nmodels = 1000
        nchains = range(6,7) 
        p_ress = range(5,6)
        p_resind = self.new_index()
        counter = 0
        #init random number generator to get reproducible random phases
        rng = random.Random(123456)
        for p_res in p_ress:
            p['p_res'] = p_res
            nchainind = self.new_index()
            for nchain in nchains:
                p['nchain'] = nchain
                #make a bunch of realizations of the same model, with randomized phases
                seqind = self.new_index()
                for i in range(0,nmodels):
                    p['seq'] = i
                    # set up chain in equal mass
                    a = [1.0]
                    for i in range(1,p['nchain']):
                        a.append(a[i-1]*((p['p_res']+1.0)/p['p_res'])**(2./3.)*rng.uniform(1.01,1.5))
                    p['a'] = a
                    p['tdep'] = 10**(rng.uniform(3,6))
                    targettime = p['tdep'] + p['deltatdep'] + 1e4*keplertime
                    #print('{} p_res {}, nchain {} seq {}'.format(counter, p['p_res'], p['nchain'], p['seq']))
                    c = chaincalc.ATauDampModel(p)
                    self.add_model(c, targettime)
                    seqind[i] = c.hash
                    counter += 1
                nchainind[nchain] = {'seq':seqind}
            p_resind[p_res] = {'nchain':nchainind}
        self.indexes['p_res'] = p_resind
