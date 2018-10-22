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

L_init_euler = data1[2][0]
E_init_euler = data1[3][0]
L_init_verlet = data2[2][0]
E_init_verlet = data2[3][0]

L_euler = data1[2][1:]
E_euler = data1[3][1:]
L_verlet = data2[2][1:]
E_verlet = data2[3][1:]

max_dev_L_euler_ind = np.argmax(abs(L_init_euler - L_euler))
max_dev_E_euler_ind = np.argmax(abs(E_init_euler - E_euler))
max_dev_L_verlet_ind = np.argmax(abs(L_init_verlet - L_verlet))
max_dev_E_verlet_ind = np.argmax(abs(E_init_verlet - E_verlet))

max_dev_L_euler = L_euler[max_dev_L_euler_ind]
max_dev_E_euler = E_euler[max_dev_E_euler_ind]
max_dev_L_verlet = L_verlet[max_dev_L_verlet_ind]
max_dev_E_verlet = E_verlet[max_dev_E_verlet_ind]


F = 15  #Fontsize
N = len(x_euler)

"""Analytical solution t values"""
t_analytical = np.linspace(0, 2*np.pi, N)
x_analytical = np.cos(t_analytical)
y_analytical = np.sin(t_analytical)

r_analytical = np.zeros(N-1)
r_euler = np.zeros_like(r_analytical)
r_verlet = np.zeros_like(r_analytical)

for i in range(N-1):
    r_analytical[i] = np.linalg.norm( [ x_analytical[i+1], y_analytical[i+1] ] )
    r_euler[i] = np.linalg.norm( [x_euler[i+1], y_euler[i+1]] )
    r_verlet[i] = np.linalg.norm( [x_verlet[i+1], y_verlet[i+1]] )

max_abs_error_euler = max((abs(r_analytical - r_euler))/r_analytical)
max_abs_error_verlet = max(abs((r_analytical - r_verlet))/r_analytical)

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

print "---------------------------------------------------"
print "N = %d  Max abs error Euler = %g Max abs error Verlet = %g " % (N, max_abs_error_euler, max_abs_error_verlet)
print "---------------------------------------------------"
print "Max_dev L     Euler = %g from init %g" % (max_dev_L_euler, L_init_euler)
print "Max_dev E     Euler = %g from init %g" % (max_dev_E_euler, E_init_euler)
print "Max_dev L    Verlet = %g from init %g" % (max_dev_L_verlet, L_init_verlet)
print "Max_dev E    Verlet = %g from init %g" % (max_dev_E_verlet, E_init_verlet)
