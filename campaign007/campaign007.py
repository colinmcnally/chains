# Colin McNally 2019 <colin@colinmcnally.ca>
# Campaign007 definitions
#  Need several tweaks to make Mp=3.16227766e-04 work
# 
import numpy as np
import sys
import chaincalc

keplertime = 2.0*np.pi*0.1**1.5
targettime = 1e9*keplertime
wall_check_interval = 15*60

class Campaign007(chaincalc.CampaignBase):
  def __init__(self):
      p = {'aspectratio0':0.035,
         'sigma0':3.8e-4*0.5*1.4,
         'redge':0.1,
         'a0':0.102,
         'deltaredge':0.001,
         'alpha':1.5,
         'flaringindex':2.0/7.0,
         'ffudge':1.0/100.0,
         'tdep':2e6*keplertime,
         'deltatdep':2e6*keplertime,
         'pmass':3.16227766e-05,
         'nchain':4,
         'q_res':1,
         'p_res':4,
         'seq':0,
         'collision':'direct',
         'G':1.0,
         'starmass':1.0,
         'integrator':'WHFAST',
         'integrator_dt':1e-2*keplertime,
         'snap_wall_interval':15*60,
         'incscatter':np.pi/100.0,
         'aspread':0.01 }
      p['physical_outputs'] = [p['tdep'], p['tdep']+p['deltatdep']]
      # try an empirical scaling of aspread
      p['aspread'] = (0.04)*(np.log(p['pmass'])-np.log(1e-6))/(np.log(1e-5) -np.log(1e-6)) + 0.01
      nmodels = 10
      nchains = range(5,11) 
      p_ress = range(3,5)
      p_resind = self.new_index()
      counter = 0
      for p_res in p_ress:
          p['p_res'] = p_res
          nchainind = self.new_index()
          for nchain in nchains:
              p['nchain'] = nchain
              #make a bunch of realizations of the same model, with randomized phases
              seqind = self.new_index()
              for i in range(0,nmodels):
                  p['seq'] = i
                  print('{} p_res {}, nchain {} seq {}'.format(counter, p['p_res'], p['nchain'], p['seq']))
                  c = chaincalc.Model(p)
                  self.add_model(c)
                  seqind[i] = c.hash
                  counter += 1
              nchainind[nchain] = {'seq':seqind}
          p_resind[p_res] = {'nchain':nchainind}
      self.indexes['p_res'] = p_resind

