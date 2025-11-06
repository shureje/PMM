import numpy as np
import matplotlib.pyplot as plt
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
plt.show()