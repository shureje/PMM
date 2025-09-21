import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Чтение данных (твой код остается)
data = []
with open('wave_data.txt', 'r') as f:
    lines = f.readlines()


current_time = None
x_vals = []
u_vals = []

for line in lines:
    line = line.strip()
    if line.startswith('t='):
        if current_time is not None:
            data.append((current_time, x_vals.copy(), u_vals.copy()))
        current_time = float(line.split('=')[1])
        x_vals.clear()
        u_vals.clear()
    elif line and not line.startswith('t='):
        parts = line.split()
        if len(parts) == 2:
            x_vals.append(float(parts[0]))
            u_vals.append(float(parts[1]))

if current_time is not None:
    data.append((current_time, x_vals.copy(), u_vals.copy()))

# Анимация
fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot([], [], 'b-', linewidth=2)
ax.set_xlim(0, 1)
ax.set_ylim(-0.5, 0.5)
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')
ax.grid(True)

def animate(frame):
    t, x, u = data[frame % len(data)]
    line.set_data(x, u)
    ax.set_title(f'Волновое уравнение, t = {t:.3f}')
    return line,

ani = animation.FuncAnimation(fig, animate, frames=len(data), 
                             interval=100, repeat=True)
plt.show()
