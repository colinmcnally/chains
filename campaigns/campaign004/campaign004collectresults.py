# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign004 collect results from outputs
# - might as well abstract this
# - dump everything from the status files to a Pandas dataframe
# 
import sys
import time
import json
import pandas as pd
import numpy as np
import chaincalc



from campaign004 import *

#Now have name col from campaign004
col = Campaign004()
chaincalc.collect_to_lifetable(col, '../lifetimes.pdh5', 'campaign004')

