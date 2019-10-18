# Colin McNally <colinmcnally.ca>
# Base class for campaigns of simulations, used by drivers and analysis to programmatically 
# define the data set
#
from collections import OrderedDict

class CampaignBase:
   """Base class for sets of models, called a Campaign
      Holds lists of models objects for now, indexed by hash
      Defines how to concatenate together, which might 
      be useful when analysing sets of results together"""
   # hold the hash-indexed Model objects
   _models = OrderedDict()
   _keys = []
   # holds dicts mapping parameter values to hases, which can then be use to lookup
   indexes = OrderedDict()

   def __add__(self, other):
     print('not implemented yet')

   def add_model(self, model):
       self._models[model.hash] = model
       self._keys.append(model.hash)       

   def new_index(self):
       return OrderedDict()

   def get_model_hash(self, key):
       return self._models[key]

   def get_model(self, index):
       return self._models[self._keys[index]]

