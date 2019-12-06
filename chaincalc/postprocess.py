# Colin McNally 2019 <colin@colinmcnally.ca>
#
# collect results from outputs of a campaign
# convert simulationarchive to recarry for orbits, and check resonance
#
import json
import numpy as np
import pandas as pd


def bool_nan_to_threestates(df, names):
    for name in names:
        df[name] = df[name].fillna(2)
        df[name] = pd.to_numeric(df[name])
    return df


def collect_to_lifetable(col, lifetable_filename, group_name):
    print('Campaign has {} models'.format(col.get_size()))

    # Don't handle this array right currently. Just cut it out.
    excludekeys = ['physical_outputs']

    data_rows = {}
    for runi in range(0, col.get_size()):
        tm = col.get_model(runi)
        targettime = col.get_targettime(runi)
        try:
            with open(tm.status_filename,'r') as json_file:
                status = json.load(json_file)
                row = {'hash': tm.hash}
                row['targettime'] = targettime
                for key in status:
                    if key=='params':
                        for pkey in status['params']:
                            if pkey in excludekeys:
                                pass
                            elif type(status['params'][pkey]) == list:
                                for i,v in enumerate(status['params'][pkey]):
                                    row[str(pkey)+str(i)] = v
                            else:
                                row[pkey] = status['params'][pkey]
                    elif type(status[key]) == list:
                        #print(key,status[key],type(status[key])==list)
                        for i,v in enumerate(status[key]):
                            #print(v,type(v))
                            row[str(key)+str(i)] = v
                    else:
                        row[key] = status[key]
                data_rows[tm.hash] = row               

        except FileNotFoundError as fe:
            print('File not found ',fe.args[0])


    lifetimes = pd.DataFrame.from_dict(data_rows, orient='index')
    strs = ['hash','status','collision','integrator','collided']
    # replace mix bool and NaN columns with 0,1,2 int 
    #lifetimes['is_this_a_constant_ratio_chain0'] = lifetimes['is_this_a_constant_ratio_chain0'].fillna(2)
    #lifetimes['is_this_a_constant_ratio_chain0'] = pd.to_numeric(lifetimes['is_this_a_constant_ratio_chain0'])
    lifetimes = bool_nan_to_threestates(lifetimes, ['is_this_a_constant_ratio_chain0', 'is_this_a_constant_ratio_chain1'])
    print(lifetimes['is_this_a_constant_ratio_chain0'].head())
    #lifetimes = lifetimes.rename({'is_this_a_constant_ratio_chain0':'is_this_a_constant_ratio_chain0_orig'})
    for i in lifetimes.dtypes.index:
        if lifetimes[i].dtype == object:
            lifetimes[i] = lifetimes[i].fillna('')
    print(lifetimes.head())

    lifetimes.to_hdf(lifetable_filename, group_name, mode='a', format='table')
    return lifetimes


def wrap2pi(x):
    return np.mod(x + np.pi, 2.0*np.pi) - np.pi


class OrbitArray:
    def __init__(self, sa):
        """ Take a rebound simulationarchive and extract a big array of orbit data """
        self.nchain = len(sa[0].particles) -1

        orbitdata = np.zeros([self.nchain,len(sa)], dtype=[('t',float), ('a',float), ('P',float), ('l',float),
            ('Omega',float), ('omega',float), ('pomega',float), ('e',float)])
        for i, sim in enumerate(sa):
            orbits = sim.calculate_orbits() # Jacobi coordinates, for heliocentraic do primary=sim.particles[0]
            for ip in range(0,self.nchain):
                orbitdata['t'][ip,i] = sim.t
                orbitdata['a'][ip,i] = sim.particles[ip+1].a
                orbitdata['P'][ip,i] = sim.particles[ip+1].P
                orbitdata['l'][ip,i] = sim.particles[ip+1].l
                orbitdata['omega'][ip,i] = sim.particles[ip+1].omega
                orbitdata['Omega'][ip,i] = sim.particles[ip+1].Omega
                orbitdata['pomega'][ip,i] = orbits[ip].pomega
                orbitdata['e'][ip,i] = orbits[ip].e
        self.orbitdata = orbitdata


    def tail_orbitarray(self, tlook):
        """Like tail, but for the orbit data arrays"""
        ibegin = np.argmax(self.orbitdata['t'][0,:] >= self.orbitdata['t'][0,-1]-tlook)
        return self.orbitdata[:][:,ibegin:]


    def compresangles(self, inner, outer, p):
        """First order resonant angles, both references"""
        #varphi = Omega + omega
        outp = p+1
        theta1 = wrap2pi( outp*outer['l'] -p*inner['l']
                      - (outer['Omega'] +outer['omega']) )
        theta2 = wrap2pi( outp*outer['l'] -p*inner['l']
                      - (inner['Omega'] +inner['omega']) )
        return [theta1, theta2]


    def compute_tight_angles(self, lookorbit, librationcut=0.90*2.0*np.pi, pmin=1, pmax=8, pratiocut=0.1):
        """Look for librating first order resonant angles, with cut off range librationcut and period ratio
           with relative cut at pratiocut"""
        tightangles = []
        for ip in range(0,self.nchain-1):
            ptight = []
            for p in range(pmin,pmax):
                angles = self.compresangles(lookorbit[:][ip,:], lookorbit[:][ip+1,:], p)
                for angle in angles:
                    angle = wrap2pi(angle-angle.mean())
                    #print(ip,ip+1, p,'range',angle.max()-angle.min(), librationcut )
                    anglerange = angle.max()-angle.min()
                    if anglerange < librationcut:
                        pratio = (lookorbit['P'][ip+1,:]/lookorbit['P'][ip,:]).mean()
                        if ((p+1.0)/p)*(1.0+pratiocut) > pratio and ((p+1.0)/p)*(1.0-pratiocut) < pratio:
                            ptight.append({'p_res':p,'rangle':anglerange, 'Pratio':pratio})
            tightangles.append(ptight)
        return tightangles


    def is_this_a_constant_ratio_chain(self, tightangles):
        """Determine from the output of compute_tight_angles if this is a constant period ratio resonant chain"""
        for ip in range(0,self.nchain-1):
            if len(tightangles[ip]) > 0:
                p = tightangles[ip][0]['p_res']
                rangle = tightangles[ip][0]['rangle']
                for ia in range(1,len(tightangles[ip])):
                    if tightangles[ip][ia]['rangle'] < rangle:
                        p = tightangles[ip][ia]['p_res']
                        rangle =  tightangles[ip][ia]['rangle']
                if ip==0:
                    masterp = p
                else:
                   if p != masterp:
                       print("ip ",ip,"p",p,"masterp",masterp)
                       return (False, None)
            else:
                print("no tight angles for ",ip)
                return (False, None)
        return (True, masterp)

