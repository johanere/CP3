import numpy as np
"""Program computes maximum absolute error of total angular
momentum for M=3 bodies, NB: filename must be modified for
specific N,T"""


M = 9   #Number of bodies
textfile = np.loadtxt("9planets_100000stepspryr_1yrs.txt")
data = textfile[1::, ::]

L1 = np.zeros(len(data[:, 0]) - M)
L2 = np.zeros_like(L1)
L3 = np.zeros_like(L1)
L4 = np.zeros_like(L1)
L5 = np.zeros_like(L1)
L6 = np.zeros_like(L1)
L7 = np.zeros_like(L1)
L8 = np.zeros_like(L1)
L9 = np.zeros_like(L1)


for i in range(0, len(data[:,0])-M):
    row = data[i, :]
    L1[i]   = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+1, :]
    L2[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+2, :]
    L3[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+3, :]
    L4[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+4, :]
    L5[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+5, :]
    L6[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+6, :]
    L7[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+7, :]
    L8[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+8, :]
    L9[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)


L_tot = L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9

L_tot_init = L1[0] + L2[0] + L3[0] + L4[0] + L5[0] + L6[0] + \
L7[0] + L8[0] + L9[0]



abs_err_L_tot = np.max(abs(L_tot - L_tot_init)/abs(L_tot_init))
print abs_err_L_tot
