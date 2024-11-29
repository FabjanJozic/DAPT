import numpy as np

g = []
r = []

data_file = 'RDF_rho0.8_T1.0_N108'
rho = 0.8
T = 1.0
N = 108

with open(data_file, 'r') as data:
    lin = data.readlines()
    for l in range(len(lin)):
        rval, gval = lin[l].strip().split()
        r.append(rval)
        g.append(gval)
    data.close()

integ = 0.0
dr = r[1]-r[0]
for i in range(len(r)):
    U_LJ = 4*(r[i]**(-12)-r[i]**(-6))
    integ += g[i]*U_LJ*(r[i]**2)*dr

U = 2*np.pi*rho*N*integ+1.5*N*T

print(U)