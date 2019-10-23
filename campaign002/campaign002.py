# Colin McNally 2019 <colin@colinmcnally.ca>
# Campaign002 definitions
# 
import numpy as np
import sys
import chaincalc

keplertime = 2.0*np.pi*0.1**1.5
targettime = 1e9*keplertime
wall_check_interval = 15*60


class Campaign002(chaincalc.CampaignBase):
  def __init__(self):
      p = {
         'tau_a':-1e6*keplertime,
         'tau_e':-1e3*keplertime,
         'tau_inc':-1e3*keplertime,
         'redge':0.1,
         'deltaredge':1.0/100.0,
         'tdep':5e5*keplertime,
         'deltatdep':5e5*keplertime,
         'a0':0.11,
         'pmass':1e-5,
         'nchain':4,
         'q_res':1,
         'p_res':4,
         'seq':0,
         'collision':'line',
         'G':1.0,
         'starmass':1.0,
         'integrator':'WHFAST',
         'integrator_dt':1e-2*keplertime,
         'snap_wall_interval':15*60,
         'incscatter':np.pi/100.0,
         'aspread':0.05 }
      p['physical_outputs'] = [p['tdep'], p['tdep']+p['deltatdep']]
      nmodels = 10
      nchains = range(5,11) 
      p_ress = range(3,7)
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
                  c = chaincalc.TauDampModel(p)
                  self.add_model(c)
                  seqind[i] = c.hash
                  counter += 1
              nchainind[nchain] = {'seq':seqind}
          p_resind[p_res] = {'nchain':nchainind}
      self.indexes['p_res'] = p_resind

