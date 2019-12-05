# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign009 collect results from outputs
# - dump everything from the status files to a Pandas dataframe
# 
import chaincalc
from campaign009 import *
col = Campaign009()
table = chaincalc.collect_to_lifetable(col, './statustable_campaign009.pdh5', 'campaign009')
