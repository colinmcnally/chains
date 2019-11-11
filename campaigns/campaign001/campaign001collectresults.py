# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign001 collect results from outputs
# - might as well abstract this
# - dump everything from the status files to a Pandas dataframe
# 
import sys
import time
import json
import pandas as pd
import numpy as np
import chaincalc


from campaign001 import *

col = Campaign001()
print('Campaign has {} models'.format(col.get_size()))
chaincalc.collect_to_lifetable(col, '../lifetimes.pdh5', 'campaign001')
