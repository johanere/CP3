import numpy as np
import matplotlib.pyplot as plt

filename1 = "EarthSunPositionsEulerMethod.dat"
filename2 = "EarthSunPositionsVerletMethod.dat"

data1 = np.loadtxt(filename1, unpack=True)
data2 = np.loadtxt(filename2, unpack=True)
x_euler = data1[0]
y_euler = data1[1]
x_verlet = data2[0]
y_verlet = data2[1]

F = 15  #Fontsize
N = len(x_euler)

"""Analytical solution t values"""
t_analytical = np.linspace(0, 2*np.pi, N)
x_analytical = np.cos(t_analytical)
y_analytical = np.sin(t_analytical)


plt.figure(figsize=(8,8))
plt.title("Two-Body system, Earth orbit around the sun, N=%d, step size %g" % (N, 1./N))
plt.plot(x_euler, y_euler, label="Forward Euler Algorithm")
plt.plot(x_verlet, y_verlet, label="Velocity-Verlet Algorithm")
plt.plot(x_analytical, y_analytical, "--", label="Analytical Solution")
plt.plot(0, 0, "ro", label="Sun")
plt.plot(x_euler[0], y_euler[0], "bo", label="Earth initial position")
plt.xlabel("x", fontsize=F)
plt.ylabel("y", fontsize=F)
plt.grid(True)
plt.legend(loc="best")
plt.show()
