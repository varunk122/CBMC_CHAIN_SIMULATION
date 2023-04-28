import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

#import all parameters
from config import *

def visualize(positions):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for pos in positions:
        ax.plot([pos[0][0],pos[1][0]],[pos[0][1],pos[1][1]],
        [pos[0][2],pos[1][2]], 'g')
        ax.scatter3D([pos[0][0]],[pos[0][1]], [pos[0][2]], color = 'red')
        ax.scatter3D([pos[1][0]],[pos[1][1]], [pos[1][2]], color = 'blue')
    plt.show()

def visualize_propane(positions):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for pos in positions:
        ax.plot([pos[0][0],pos[1][0],pos[2][0]],[pos[0][1],pos[1][1],pos[2][1]],
        [pos[0][2],pos[1][2],pos[2][2]], 'g')
        ax.scatter3D([pos[0][0]],[pos[0][1]], [pos[0][2]], color = 'red')
        ax.scatter3D([pos[1][0]],[pos[1][1]], [pos[1][2]], color = 'blue')
        ax.scatter3D([pos[2][0]],[pos[2][1]], [pos[2][2]], color = 'yellow')
    plt.show()

def visualize_butane(positions):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for pos in positions:
        ax.plot([pos[0][0],pos[1][0],pos[2][0],pos[3][0]],[pos[0][1],pos[1][1],pos[2][1],pos[3][1]],
        [pos[0][2],pos[1][2],pos[2][2],pos[3][2]], 'g')
        ax.scatter3D([pos[0][0]],[pos[0][1]], [pos[0][2]], color = 'red')
        ax.scatter3D([pos[1][0]],[pos[1][1]], [pos[1][2]], color = 'blue')
        ax.scatter3D([pos[2][0]],[pos[2][1]], [pos[2][2]], color = 'yellow')
        ax.scatter3D([pos[3][0]],[pos[3][1]], [pos[3][2]], color = 'purple')
    plt.show()

def plot_energy_from_file(file_name):
    with open(file_name,'r') as f:
        energy = [float(line.strip()) for line in f]
        steps = [ int(write_interval * i) for i in range(len(energy))]
        plt.plot(steps, energy)
        plt.grid()
        plt.xlabel("Steps")
        plt.ylabel("Energy (KbT)")
        plt.show()

def convert_positions_from_file_to_video(positions_list):
    
    def animate(i):
        positions = positions_list[i]
        # print(positions_list)
        for pos in positions:
            # ax.plot([pos[0][0],pos[1][0]],[pos[0][1],pos[1][1]],
            # [pos[0][2],pos[1][2]], 'g')
            # ax.scatter3D([pos[0][0]],[pos[0][1]], [pos[0][2]], color = 'red')
            # ax.scatter3D([pos[1][0]],[pos[1][1]], [pos[1][2]], color = 'blue')
            ax.plot([pos[0][0],pos[1][0],pos[2][0],pos[3][0]],[pos[0][1],pos[1][1],pos[2][1],pos[3][1]],
            [pos[0][2],pos[1][2],pos[2][2],pos[3][2]], 'g')
            ax.scatter3D([pos[0][0]],[pos[0][1]], [pos[0][2]], color = 'red')
            ax.scatter3D([pos[1][0]],[pos[1][1]], [pos[1][2]], color = 'blue')
            ax.scatter3D([pos[2][0]],[pos[2][1]], [pos[2][2]], color = 'yellow')
            ax.scatter3D([pos[3][0]],[pos[3][1]], [pos[3][2]], color = 'purple')

    num_plots = len(positions_list)
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection='3d')
    ani = FuncAnimation(fig, animate, frames=num_plots, interval=200)
    writer = PillowWriter(fps=4)
    ani.save('plots.gif', writer=writer)

