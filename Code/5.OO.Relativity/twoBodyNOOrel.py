import numpy as np
import os
import matplotlib.pyplot as plt

data = np.loadtxt("MercurySunPositionsVerletMethod.dat", unpack=True)
x = data[0]
y = data[1]
z = data[2]
j_=5
T = 100
F = 13



plt.figure(figsize=(8,8))
plt.title(r"Two-Body system, Mercury-Sun, $n=10^%d$, T=%d yrs" % (j_, T ))
plt.plot(x, y, label="Mercury$")
plt.plot(0, 0, "ro", label="Sun")
plt.xlabel("x [AU]", fontsize=F)
plt.ylabel("y [AU]", fontsize=F)
plt.axis("equal")
plt.tick_params(axis='both', which='major', labelsize=F)
plt.grid(True)
plt.legend(loc="upper right", fontsize=F-1)
plt.show()

Mercury_perihelion_distance = 0.3075
tol = 1.e-6
x_p = x[-1]
y_p = y[-1]
conv_to_arcsec = 3600.0
perihelions = []
#Find last perihelion
for i in range(len(x)):
    dist = np.sqrt(x[i]**2 + y[i]**2)
    if abs(dist - Mercury_perihelion_distance) < tol:
        perihelions.append(i)
    else:
        pass

last = perihelions[-1]
x_p = x[last]
y_p = y[last]
periangle = np.arctan(abs(y_p/x_p))*conv_to_arcsec
print periangle
