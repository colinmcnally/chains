# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign007 collect results from outputs
# - dump everything from the status files to a Pandas dataframe
# 
import chaincalc
from campaign007 import *
col = Campaign007()
chaincalc.collect_to_lifetable(col, '../lifetimes.pdh5', 'campaign007')

