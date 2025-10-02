import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def load_all_data():
    data = []
    filepath = os.path.join(os.path.dirname(__file__), 'transport_data.dat')
    with open(filepath, 'r') as f:
        current_data = []
        for line in f:
            if line.startswith('# Step'):
                if current_data:
                    arr = np.array(current_data)
                    x = arr[:, 0].reshape(100, 100)
                    y = arr[:, 1].reshape(100, 100)
                    rho = arr[:, 2].reshape(100, 100)
                    data.append((x, y, rho))
                current_data = []
            elif line.strip() and not line.startswith('#'):
                vals = line.strip().split()
                if len(vals) == 3:
                    current_data.append([float(v) for v in vals])
        
        if current_data:
            arr = np.array(current_data)
            x = arr[:, 0].reshape(100, 100)
            y = arr[:, 1].reshape(100, 100)
            rho = arr[:, 2].reshape(100, 100)
            data.append((x, y, rho))
    
    return data

def create_animation():
    data = load_all_data()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    def animate(frame):
        ax.clear()
        x, y, rho = data[frame]
        
        print(f"Frame {frame}: min={np.min(rho):.6f}, max={np.max(rho):.6f}")
        
        im = ax.contourf(x, y, rho, levels=50, cmap='viridis')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Перенос примесей в 2D (шаг {frame})')
        return [im]
    
    anim = FuncAnimation(fig, animate, frames=len(data), interval=50, repeat=True)
    plt.show()

if __name__ == "__main__":
    create_animation()
