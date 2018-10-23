
import numpy as np
import matplotlib.pyplot as plt
import string as str
outdata=["2_data_10.txt"]


for j in range(len(outdata)):  
    data=np.loadtxt(outdata[j],unpack=True)
    paramters=data[:,0]
    N=int(paramters[0])
    printstep=paramters[1]
    M=int(paramters[2]) #numb of planets
    
    L=np.zeros((N,3))

    Ltemp=np.zeros(3)

    for j in range(1,M*N,M):
        for i in range(0,M):
            Ltemp+=data[7:,i+j]  
        L[i+j/2-1]=Ltemp
        Ltemp[:]=np.zeros(3)
 

    Ltot=np.sum(L, axis=1)
    
    print N
    time=np.linspace(0,1,N)

    print len(time)
    print len(Ltot)
    
#===============================================================================
#     for i in range(0,M):
#         plt.plot(data[1,1+i::M],data[2,1+i::M],label='%s'%planets[i])
#     
# plt.xlabel('$x[AU]$')
# plt.ylabel('$y[AU]$')
# plt.grid('on')
# plt.title('The Solar system with %s bodies'%M)
# plt.legend(loc=0)
# plt.show()       
