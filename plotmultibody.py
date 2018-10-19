import numpy as np
import matplotlib.pyplot as plt
import string as str
outdata=["3_data_1000.txt"]

planets=["Sun","Earth","Mars"]

for j in range(len(outdata)):  
    data=np.loadtxt(outdata[j],unpack=True)
    n=int(data[0,0])
    m=int(data[0,1])
    
    pos=np.zeros((n,m))
    vel=np.zeros((n,m))


    pos=np.asarray(data[:,2::2])
    vel=np.asarray(data[:,3::2])
    
    xpos=pos[0]
    ypos=pos[1]
    zpos=pos[2]

    xvel=vel[0]
    yvel=vel[1]
    zvel=vel[2]
 
    for i in range(0,m):
        plt.plot(xpos[i::m],ypos[i::m],label='%s'%planets[i])
         
plt.xlabel('$x[AU]$')
plt.ylabel('$y[AU]$')
plt.grid('on')
plt.title('The Solar system with %s bodies'%m)
plt.legend(loc=0)
plt.show()
