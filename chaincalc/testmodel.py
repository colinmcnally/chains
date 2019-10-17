import chainmodel

p = {'aspectratio0':0.035,
     'sigma0':3.8e-4,
     'redge':0.1,
     'deltaredge':0.001,
     'alpha':1.0,
     'flaringindex':2.0/7.0,
     'ffudge':1.0/100.0,
     'tdep':1e5,
     'deltatdep':1e3,
     'pmass':1e-5,
     'nchain':9,
     'p_res':4 }
c = chainmodel.model(p)

print(c.hash)
