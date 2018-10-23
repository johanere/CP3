import numpy as np
import sys
import matplotlib.pyplot as plt

filename1 = "EarthSunPositionsEulerMethod.dat"
filename2 = "EarthSunPositionsVerletMethod.dat"

printstep = float(sys.argv[1])

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

max_abs_err_euler_L  = max(abs(L_init_euler - L_euler)/L_init_euler)
max_abs_err_euler_E  = max(abs(E_init_euler - E_euler)/abs(E_init_euler))
max_abs_err_verlet_L = max(abs(L_init_verlet - L_verlet)/L_init_verlet)
max_abs_err_verlet_E = max(abs(E_init_verlet - E_verlet)/abs(E_init_verlet))


F = 12  #Fontsize
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
    r_euler[i]      = np.linalg.norm( [ x_euler[i+1], y_euler[i+1] ] )
    r_verlet[i]     = np.linalg.norm( [ x_verlet[i+1], y_verlet[i+1] ] )

max_abs_error_euler_r  = max(abs(r_analytical - r_euler )/r_analytical)
max_abs_error_verlet_r = max(abs(r_analytical - r_verlet)/r_analytical)

plt.figure(figsize=(8,8))
plt.title("Two-Body system, Earth orbit around Sun fixed at origin, N=%d, step size h= %g" % (N*printstep, 1./(N*printstep)) )
plt.plot(x_euler, y_euler, label="Forward Euler Algorithm")
plt.plot(x_verlet, y_verlet, label="Velocity-Verlet Algorithm")
plt.plot(x_analytical, y_analytical, "--", label="Analytical Solution")
plt.plot(0, 0, "ro", label="Sun")
plt.plot(x_euler[0], y_euler[0], "bo", label="Earth initial position")
plt.xlabel("x", fontsize=F)
plt.ylabel("y", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.grid(True)
plt.legend(loc="upper right", fontsize=F-1)
plt.show()


print "----------------------------------------------------------"
print "|    N = %d    |      Euler        |       Verlet      |" % N
print "----------------------------------------------------------"
print "|    eps_r     |      %-8.6g     |       %-8.6g | " % (max_abs_error_euler_r, max_abs_error_verlet_r)
print "|    eps_E     |      %-8.6g   |       %-8.6g    | " % (max_abs_err_euler_E, max_abs_err_verlet_E)
print "|    eps_L     |      %-8.6g     |       %-8.6g    | " % (max_abs_err_euler_L, max_abs_err_verlet_L)
print "----------------------------------------------------------"
