import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Чтение данных
with open('cmd/result.txt', 'r') as f:
    data = [list(map(float, line.split())) for line in f if line.strip()]

# data[i][j] - i-я точка пространства, j-й момент времени
N = len(data)  # количество точек по пространству
M = len(data[0])  # количество моментов времени

fig, ax = plt.subplots()
x = np.linspace(0, 1, N)  

ax.set_xlim(0, 1)
ax.set_ylim(0, 3)
ax.set_xlabel('x')
ax.set_ylabel('T(x,t)')
ax.set_title('Уравнение теплопроводности')
ax.grid(True)

line, = ax.plot([], [], 'b-', linewidth=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def animate(frame):
    y = [data[i][frame] for i in range(N)]
    line.set_data(x, y)
    time_text.set_text(f'Время: {frame*0.005:.3f}')
    return line, time_text

anim = animation.FuncAnimation(fig, animate, frames=M, 
                             interval=50, blit=True, repeat=True)
plt.show()
