import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt('heat_animation.txt')
x = data[:, 0]
y = data[:, 1]
t = data[:, 2].astype(int)
u = data[:, 3]

nx = 50
ny = 50
time_steps = int(len(data) / (nx * ny))

u_animation = []
for step in range(time_steps):
    start = step * nx * ny
    end = start + nx * ny
    u_step = u[start:end].reshape(nx, ny)
    u_animation.append(u_step)

vmin = min(np.min(u_arr) for u_arr in u_animation)
vmax = max(np.max(u_arr) for u_arr in u_animation)

fig, ax = plt.subplots(figsize=(10, 6))

#цветовая полоса
contour = ax.contourf(u_animation[0], levels=20, extent=[0, 2, 0, 3], 
                     cmap='viridis', vmin=vmin, vmax=vmax)
cbar = fig.colorbar(contour, ax=ax)
cbar.set_label('Температура')

def animate(frame):
    for coll in ax.collections:
        coll.remove()
    
    contour = ax.contourf(u_animation[frame], levels=20, extent=[0, 2, 0, 3], 
                         cmap='viridis', vmin=vmin, vmax=vmax)
    ax.contour(u_animation[frame], levels=20, extent=[0, 2, 0, 3], 
               colors='black', linewidths=0.5)
    t = (frame * 1.5/time_steps)
    ax.set_title(f'Время: {t:.2f}')
    return contour,

anim = FuncAnimation(fig, animate, frames=time_steps, interval=100, blit=False)
plt.show()