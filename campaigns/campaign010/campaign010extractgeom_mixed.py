# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Extract chain geometry, mixed ratio chains in campaign010
# 
import pandas as pd
import rebound
import chaincalc
from campaign010 import *
import matplotlib.pyplot as plt

st = pd.read_hdf('statustable_campaign010.pdh5','campaign010')

#Some frequent chains
#ps = [3, 4, 4, 4, 4]
#ps = [3, 3, 4, 4, 4]
#ps = [3, 3, 3, 4, 4]
#ps = [3, 3, 4, 4, 3]
#ps = [4, 4, 4, 4, 3]
#ps = [4, 5, 4, 5, 4]
ps = [5, 4, 5, 4, 5]
#ps = [3, 3, 4, 5, 5]
#ps = [4, 3, 4, 5, 5]
#ps = [3, 3, 4, 5, 4]

#ps = [4, 3, 4, 3, 4]
  
set_in_target = (st.is_this_a_mixed_ratio_chain1_0_0==ps[0]) & (st.is_this_a_mixed_ratio_chain1_1_0==ps[1]) & (st.is_this_a_mixed_ratio_chain1_2_0==ps[2]) & (st.is_this_a_mixed_ratio_chain1_3_0==ps[3]) & (st.is_this_a_mixed_ratio_chain1_4_0==ps[4])  


def limit_to_ten(series, limit = 10):
    nt = 0
    for index, val in series.iteritems():
        if val:
            nt += 1
            if nt > 10:
                series[index] = False
    return series



def get_longitude_of_periastron_stats(row):
    print('Processing ',row.simarchive_filename)
    sa = rebound.SimulationArchive(row.simarchive_filename, process_warnings=False)
    oa = chaincalc.OrbitArray(sa)
    tlook = 1e4*keplertime #probably 1e2-1e3 times the output cadence
    lko = oa.tail_orbitarray(tlook)
    mean_dlongitude_of_periastron = []
    std_dlongitude_of_periastron = []
    for n in range(lko['pomega'].shape[0]-1):
        longstats = [ {'mean':chaincalc.wrap2pi(lko['pomega'][n+1,:] -lko['pomega'][n,:] +a).mean()-a, 
                       'std': chaincalc.wrap2pi(lko['pomega'][n+1,:] -lko['pomega'][n,:] +a).std()  } for a in [0.0, np.pi]]
        if longstats[0]['std'] < longstats[1]['std']:
            mean_dlongitude_of_periastron.append(longstats[0]['mean'])
            std_dlongitude_of_periastron.append(longstats[0]['std'])
        else:
            mean_dlongitude_of_periastron.append(longstats[1]['mean'])
            std_dlongitude_of_periastron.append(longstats[1]['std'])
    #print(mean_dlongitude_of_periastron)
    #print(std_dlongitude_of_periastron)
    #if max(mean_dlongitude_of_periastron) > -np.pi/2:
    #    print(row.hash)
    #    print(mean_dlongitude_of_periastron)
    return {'delta_longitude_of_periastron':mean_dlongitude_of_periastron, 'std_delta_longitude_of_periastron':std_dlongitude_of_periastron}


#set_in_target = limit_to_ten(set_in_target)


def plot_chain_peris(set_in_res):
    delta_longitude_of_periastron = st[set_in_res].apply(get_longitude_of_periastron_stats,
                                                         axis=1, result_type='expand')
    print(delta_longitude_of_periastron)
    for index, row in delta_longitude_of_periastron.iterrows():
        plt.errorbar(np.arange(len(row.delta_longitude_of_periastron)), row.delta_longitude_of_periastron, 
                     yerr = row.std_delta_longitude_of_periastron, capsize=10, linewidth=1)
    plt.ylabel(r'$\varpi_{i+1} - \varpi_i$')
    plt.xlabel('Planet index $i$')
    plt.ylim(-2*np.pi,2*np.pi)

plt.figure()
plot_chain_peris(set_in_target)
for y in  [-1,-0.5,0.0,0.5,1]:
    plt.axhline(y=y*np.pi, color='grey', linewidth=1)
plt.title(r'{:s} chains longitude of periastron differences'.format(str(ps)))

plt.show()
