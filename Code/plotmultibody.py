
import numpy as np
import matplotlib.pyplot as plt
import string as str
outdata=["3_data_100.txt"]

planets=["Sun","Earth","Mars"]

for j in range(len(outdata)):  
    data=np.loadtxt(outdata[j],unpack=True)
    paramters=data[:,0]
    N=int(paramters[0])
    printstep=paramters[1]
    M=int(paramters[2]) #numb of planets
    
    for i in range(0,M):
        plt.plot(data[1,1+i::M],data[2,1+i::M],label='%s'%planets[i])
    
plt.xlabel('$x[AU]$')
plt.ylabel('$y[AU]$')
plt.grid('on')
plt.title('The Solar system with %s bodies'%M)
plt.legend(loc=0)
plt.show()       
   