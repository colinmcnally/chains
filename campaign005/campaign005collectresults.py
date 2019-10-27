# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign005 collect results from outputs
# - dump everything from the status files to a Pandas dataframe
# 
import chaincalc
from campaign005 import *
col = Campaign005()
chaincalc.collect_to_lifetable(col, '../lifetimes.pdh5', 'campaign005')

