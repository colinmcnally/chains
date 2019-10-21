# Colin McNally 2019 <colin@colinmcnally.ca>
# Script to test the CampaignBase usage.
# Set to run in laptop, not container
# 
import numpy as np
import sys
sys.path.append('/Users/colinm/vc/rebound')
sys.path.append('../crbx')
sys.path.append('../')
import rebound
import chaincalc

# the test driver uses these
keplertime = 2.0*np.pi*0.1**1.5
targettime = 40000*keplertime
wall_check_interval = 15

class myruns(chaincalc.CampaignBase):
  def __init__(self):
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
         'q_res':1,
         'p_res':4,
         'seq':0,
         'collision':'line',
         'G':1.0,
         'starmass':1.0,
         'integrator':'WHFAST',
         'integrator_dt':1e-2*2*np.pi*0.1**1.5,
         'snap_wall_interval':60*60,
         'incscatter':np.pi/100.0,
         'aspread':0.05  }
      #make a bunch of realizations of the same model, with randomized phases
      nmodels = 10
      nchains = range(5,11) 
      nchainind = self.new_index()
      for nchain in nchains:
          seqind = self.new_index()
          for i in range(0,nmodels):
              p['seq'] = i
              c = chaincalc.Model(p)
              print(c.hash)
              self.add_model(c)
              seqind[i] = c.hash
          nchainind[nchain] = seqind
      self.indexes['nchain'] = nchainind


print('generate a campaign')
col = myruns()
print('')
for seqs in col.indexes['nchain']:
    print('nchain',seqs)
    for i in col.indexes['nchain'][seqs]:
        print('realization ',i,' ',col.indexes['nchain'][seqs][i])

print('')
print('test get_model method')
print(col.get_model(2))

