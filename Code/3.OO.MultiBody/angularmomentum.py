import numpy as np

M = 3   #Number of bodies
textfile = np.loadtxt("3planets_1000000stepspryr_12yrs.txt")
data = textfile[1::, ::]

L1 = np.zeros(len(data[:, 0]) - M)
L2 = np.zeros_like(L1)
L3 = np.zeros_like(L1)

for i in range(0, len(data[:,0])-M):
    row = data[i, :]
    L1[i]   = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+1, :]
    L2[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

    row = data[i+2, :]
    L3[i] = np.sqrt(row[-3]**2 + row[-2]**2 + row[-1]**2)

L_tot = L1 + L2 + L3

L_tot_init = L1[0] + L2[0] + L3[0]


abs_err_L_tot = np.max(abs(L_tot - L_tot_init)/abs(L_tot_init))
print abs_err_L_tot
