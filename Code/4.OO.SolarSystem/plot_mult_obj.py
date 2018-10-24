
import numpy as np
import matplotlib.pyplot as plt
import string as str
T = 1
j_ = 5
n = 10**j_

outdata=["9planets_100000stepspryr_1yrs.txt"]

planets=["Sun","Earth","Jupiter", "Mars", "Venus", "Saturn", "Mercury", "Uranus", \
"Neptune"]

for j in range(len(outdata)):
    data=np.loadtxt(outdata[j],unpack=True)
    paramters=data[:,0]
    N=int(paramters[0])
    printstep=paramters[1]
    M=int(paramters[2]) #numb of planets
    T=int(paramters[3])

    for i in range(0,M):
        plt.plot(data[1,1+i::M],data[2,1+i::M],label='%s'%planets[i])

plt.xlabel('$x[AU]$')
plt.ylabel('$y[AU]$')
plt.grid('on')
plt.title(r'The Solar system with %s bodies. $n=10^%d$, T=%d yrs'% (M, j_, T))
plt.legend(loc=0)
plt.show()
