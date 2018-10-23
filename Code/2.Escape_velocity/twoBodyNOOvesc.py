import numpy as np
import sys
import matplotlib.pyplot as plt

filename2 = "EarthSunPositionsVerletMethod.dat"

printstep = float(sys.argv[1])

data2 = np.loadtxt(filename2, unpack=True)
x_verlet = data2[0]
y_verlet = data2[1]


F = 12  #Fontsize
N = len(x_verlet)

plt.figure(figsize=(8,8))
plt.title("Two-Body system, Earth orbit around Sun fixed at origin, N=%d, step size h= %g" % (N*printstep, 1./(N*printstep)) )
plt.plot(x_verlet, y_verlet, label="Velocity-Verlet Algorithm")
plt.plot(0, 0, "ro", label="Sun")
plt.plot(x_verlet[0], y_verlet[0], "bo", label="Earth initial position")
plt.xlabel("x", fontsize=F)
plt.ylabel("y", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.grid(True)
plt.legend(loc="upper right", fontsize=F-1)
plt.show()
