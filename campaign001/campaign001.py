# Colin McNally 2019 <colin@colinmcnally.ca>
# Campaign001 definitions
# 
import numpy as np
import sys
import chaincalc

# the test driver uses these
keplertime = 2.0*np.pi*0.1**1.5
targettime = 1e9*keplertime
wall_check_interval = 30*60*60

class Campaign001(chaincalc.CampaignBase):
  def __init__(self):
      p = {'aspectratio0':0.035,
         'sigma0':3.8e-4*0.5*1.4,
         'redge':0.1,
         'a0':0.11,
         'deltaredge':0.001,
         'alpha':1.0,
         'flaringindex':2.0/7.0,
         'ffudge':1.0/100.0,
         'tdep':4e5*keplertime,
         'deltatdep':4e5*keplertime,
         'pmass':1e-5,
         'nchain':4,
         'q_res':1,
         'p_res':4,
         'seq':0 }
      nmodels = 10
      nchains = range(5,11) 
      p_ress = range(3,7)
      p_resind = self.new_index()
      for p_res in p_ress:
          p['p_res'] = p_res
          nchainind = self.new_index()
          for nchain in nchains:
              p['nchain'] = nchain
              #make a bunch of realizations of the same model, with randomized phases
              seqind = self.new_index()
              for i in range(0,nmodels):
                  p['seq'] = i
                  c = chaincalc.Model(p)
                  self.add_model(c)
                  seqind[i] = c.hash
              nchainind[nchain] = {'seq':seqind}
          p_resind[p_res] = {'nchain':nchainind}
      self.indexes['p_res'] = p_resind

