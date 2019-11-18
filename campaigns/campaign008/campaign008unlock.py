# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Chains campaign008 reset status locks to unlocked
# 
import chaincalc
import json

from campaign008 import *
col = Campaign008()

# should move this into postprocessing or utils or something
for runi in range(0, col.get_size()):
    tm = col.get_model(runi)
    try:
        with open(tm.status_filename,'r') as json_file:
            status = json.load(json_file)
        if status['lock']:
            print('unlocking: {}'.format(tm.status_filename))
            status['lock'] = False
            with open(tm.status_filename,'w') as json_file:
                json.dump(status, json_file)
    except FileNotFoundError:
        print('No status file was found at: {}'.format(tm.status_filename))
