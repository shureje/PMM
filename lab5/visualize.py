import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('solution.dat')
x = data[:, 0]
t = data[:, 1]
u = data[:, 2]

unique_times = np.unique(t)
N = len(np.unique(x))

X = x[:N]
T = unique_times
U = u.reshape(len(unique_times), N)

fig = plt.figure(figsize=(12, 5))

T_mesh, X_mesh = np.meshgrid(T, X)
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X_mesh, T_mesh, U.T, cmap='viridis')
ax1.set_xlabel('x')
ax1.set_ylabel('t')
ax1.set_zlabel('u(x,t)')

ax2 = fig.add_subplot(122)
for i, time in enumerate(unique_times[::2]):
    idx = np.where(t == time)[0]
    ax2.plot(x[idx], u[idx], label=f't = {time:.3f}')
ax2.set_xlabel('x')
ax2.set_ylabel('u(x,t)')
ax2.legend()
ax2.grid(True)

plt.show()