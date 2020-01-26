# Colin McNally 2019 <colin@colinmcnally.ca>
#
# Extract chain geometry, just constant ratio chains
# 
import pandas as pd
import rebound
import chaincalc
from campaign010 import *
import matplotlib.pyplot as plt

st = pd.read_hdf('statustable_campaign010.pdh5','campaign010')

set_not_in_5to4 = (st.is_this_a_constant_ratio_chain0<2) & (st.is_this_a_constant_ratio_chain1!=4)
set_in_5to4 = (st.is_this_a_constant_ratio_chain0==1) & (st.is_this_a_constant_ratio_chain1==4)
set_in_3to2 = (st.is_this_a_constant_ratio_chain0==1) & (st.is_this_a_constant_ratio_chain1==2)
set_in_4to3 = (st.is_this_a_constant_ratio_chain0==1) & (st.is_this_a_constant_ratio_chain1==3)
set_in_6to5 = (st.is_this_a_constant_ratio_chain0==1) & (st.is_this_a_constant_ratio_chain1==5)

def limit_to_ten(series, limit = 10):
    nt = 0
    for index, val in series.iteritems():
        if val:
            nt += 1
            if nt > 10:
                series[index] = False
    return series


def doublecheck(row):
    sa = rebound.SimulationArchive(row.simarchive_filename, process_warnings=False)
    oa = chaincalc.OrbitArray(sa)
    tlook = 1e4*keplertime #probably 1e2-1e3 times the output cadence
    lko = oa.tail_orbitarray(tlook)
    finding = oa.is_this_a_constant_ratio_chain(oa.compute_tight_angles(lko,pratiocut=0.02))
    print(row.hash, finding)
    return finding[0]

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
    if max(mean_dlongitude_of_periastron) > -np.pi/2:
        print(row.hash)
        print(mean_dlongitude_of_periastron)
    return {'delta_longitude_of_periastron':mean_dlongitude_of_periastron, 'std_delta_longitude_of_periastron':std_dlongitude_of_periastron}


#set_in_5to4 = limit_to_ten(set_in_5to4)

#update aligns on the Series index, and overwrites values
set_in_6to5.update(st[set_in_6to5].apply(doublecheck, axis=1))
set_in_5to4.update(st[set_in_5to4].apply(doublecheck, axis=1))
set_in_4to3.update(st[set_in_4to3].apply(doublecheck, axis=1))
#set_in_3to2.update(st[set_in_3to2].apply(doublecheck, axis=1))


def plot_chain_peris(set_in_res):
    delta_longitude_of_periastron = st[set_in_res].apply(get_longitude_of_periastron_stats,
                                                         axis=1, result_type='expand')
    print(delta_longitude_of_periastron)
    for index, row in delta_longitude_of_periastron.iterrows():
        plt.errorbar(np.arange(len(row.delta_longitude_of_periastron)), row.delta_longitude_of_periastron, yerr = row.std_delta_longitude_of_periastron, capsize=10)
    plt.ylabel(r'$\varpi_{i+1} - \varpi_i$')
    plt.xlabel('Planet index $i$')
    plt.ylim(-2*np.pi,2*np.pi)

plt.figure()
plot_chain_peris(set_in_5to4)
plt.title(r'5:4 chains longitude of periastron differences')

plt.figure()
plot_chain_peris(set_in_4to3)
plt.title(r'4:3 chains longitude of periastron differences')

plt.figure()
plot_chain_peris(set_in_6to5)
plt.title(r'6:5 chains longitude of periastron differences')

#plot_chain_peris(set_in_3to2)
#plt.title(r'3:2 chains longitude of periastron differences')


plt.show()
