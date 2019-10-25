# Colin McNally 2019 <colin@colinmcnally.ca>
#
# collect results from outputs of a campaign
# 
import json
import pandas as pd


def collect_to_lifetable(col, lifetable_filename, group_name):
    print('Campaign has {} models'.format(col.get_size()))

    # Don't handle this array right currently. Just cut it out.
    excludekeys = ['physical_outputs']

    data_rows = {}
    for runi in range(0, col.get_size()):
        tm = col.get_model(runi)
        try:
            with open(tm.status_filename,'r') as json_file:
                status = json.load(json_file)
                row = {'hash': tm.hash}
                for key in status:
                    if key=='params':
                        for pkey in status['params']:
                            if pkey in excludekeys:
                                pass
                            else:
                                row[pkey] = status['params'][pkey]
                    elif type(status[key]) == list:
                        for i,v in enumerate(status[key]):
                            row[str(key)+str(i)] = v
                    else:
                        row[key] = status[key]
                data_rows[tm.hash] = row               

        except FileNotFoundError as fe:
            print('File not found ',fe.args[0])


    lifetimes = pd.DataFrame.from_dict(data_rows, orient='index')
    strs = ['hash','status','collision','integrator','collided']
    for i in lifetimes.dtypes.index:
        if lifetimes[i].dtype == object:
            lifetimes[i].fillna('')
    print(lifetimes.head())

    lifetimes.to_hdf(lifetable_filename, group_name, mode='a', format='table')
