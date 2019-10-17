from collections import OrderedDict

class CampaignBase:
   """Base class for sets of models, called a Campaign
      Holds lists of models objects for now, indexed by hash
      Defines how to concatenate together, which might 
      be useful when analysing sets of results together"""
   # hold the hash-indexed Model objects
   models = OrderedDict()
   # holds dicts mapping parameter values to hases, which can then be use to lookup
   indexes = OrderedDict

   def __add__(self, other):
     print('not implemented yet')

