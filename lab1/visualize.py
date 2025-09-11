import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Параметры
L = 1.0
N = 101
x = np.linspace(0, L, N)

# Загрузка данных
data = np.loadtxt('temperature_data.txt')

fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot(x, data[0], 'b-', linewidth=2)
ax.set_xlim(0, L)
ax.set_ylim(-0.5, 2.5)
ax.set_xlabel('x')
ax.set_ylabel('T(x,t)')
ax.set_title('Распространение температуры')
ax.grid(True)

def animate(frame):
    line.set_ydata(data[frame])
    ax.set_title(f'Распространение температуры (кадр {frame})')
    return line,

ani = animation.FuncAnimation(fig, animate, frames=len(data), interval=100, blit=False, repeat=True)
plt.show()