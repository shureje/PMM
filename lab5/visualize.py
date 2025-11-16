import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
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

=======
from matplotlib.animation import FuncAnimation

data = np.genfromtxt('solution.dat', invalid_raise=False)
data = data[~np.isnan(data).any(axis=1)]

N = 50
M = len(data) // N
total_points = M * N
data = data[:total_points]

x = data[:, 0]
t = data[:, 1]
T = data[:, 2]

X = x.reshape(M, N)
Time = t.reshape(M, N)
Temp = T.reshape(M, N)

fig, ax = plt.subplots()
line, = ax.plot([], [], 'b-', linewidth=2)
ax.set_xlim(0, 1)
ax.set_ylim(0, 4)
ax.set_xlabel('x')
ax.set_ylabel('T(x,t)')
ax.grid(True)
title = ax.set_title('')

def animate(frame):
    idx = frame * 50
    line.set_data(X[idx], Temp[idx])
    title.set_text(f'Время t = {Time[idx, 0]:.4f}')
    return line, title

ani = FuncAnimation(fig, animate, frames=M//50, interval=1, blit=True, repeat=True)
>>>>>>> 8b9d4ff873986debc7c15205c127e8e043cf5dbc
plt.show()