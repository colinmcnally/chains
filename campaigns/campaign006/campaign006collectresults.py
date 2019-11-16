# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign006 collect results from outputs
# - dump everything from the status files to a Pandas dataframe
# 
import chaincalc
from campaign006 import *
col = Campaign006()
chaincalc.collect_to_lifetable(col, '../lifetimes.pdh5', 'campaign006')

