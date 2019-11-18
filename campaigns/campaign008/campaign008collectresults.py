# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign008 collect results from outputs
# - dump everything from the status files to a Pandas dataframe
# 
import chaincalc
from campaign008 import *
col = Campaign008()
lifetable = chaincalc.collect_to_lifetable(col, './statustable.pdh5', 'campaign008')

