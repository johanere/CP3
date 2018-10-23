import matplotlib.pyplot as plt
import numpy as np
import os

M = 6                           #Number of beta values to be plotted
betas = np.linspace(2, 3, M)    #Array for M beta values in range [2, 3]
T = 1
N = 100000
printstep = 5

cppFile = "twoBodyNOOvesc.cpp"

F = 14      #Fontsize

plt.figure(figsize=(8,8))
plt.title(r"""Earth-Sun system with increasing $\beta$.
Earth initial velocity = $2\pi+1$, N=%d, T=%d, printstep=%d""" %(N, T, printstep))

for i in range(M):
    os.system("c++ -std=c++11 %s -o proj && ./proj %d %d %d %f && rm proj" % (cppFile, N, printstep, T, betas[i]))
    data = np.loadtxt("EarthSunPositionsVerletMethod.dat", unpack=True)
    plt.plot(data[0], data[1], label=r"$\beta=%.3g$" %betas[i])
plt.xlabel("x [AU]", fontsize=F)
plt.ylabel("y [AU]", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F-1)
plt.plot(0,0, "o", label="Sun")
plt.plot(1,0, "o", label="Earth Initial")
plt.legend()
plt.grid(True)
plt.show()
