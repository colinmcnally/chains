# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign010 collect results from outputs
# - dump everything from the status files to a Pandas dataframe
# 
import chaincalc
from campaign010 import *
col = Campaign010()
table = chaincalc.collect_to_lifetable(col, './statustable_campaign010.pdh5', 'campaign010')
