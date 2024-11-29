import numpy as np

data_file = 'RDF_rho0.8_T1.0_N108'
rho = 0.8
T = 1.0
N = 108

with open('C:\\Users\\fabja\\Desktop\\programi\\DAPT\\vjezba5\\'+data_file, 'r') as data:
    lin = data.readlines()
    r = [0]*len(lin)
    g = [0]*len(lin)
    dr = [0.0]
    for l in range(len(lin)):
        rval, gval = lin[l].strip().split()
        r[l] = float(rval)
        g[l] = float(gval)
        dr.append(r[l])
    data.close()

integ = 0.0
for i in range(len(r)):
    U_LJ = 4*(r[i]**(-12)-r[i]**(-6))
    integ += g[i]*U_LJ*(r[i]**2)*(dr[i+1]-dr[i])

U = 2*np.pi*rho*N*integ+1.5*N*T

print(U)