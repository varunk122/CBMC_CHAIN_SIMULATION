import matplotlib.pyplot as plt

def visualize(positions):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for pos in positions:
        ax.plot([pos[0][0],pos[1][0]],[pos[0][1],pos[1][1]],
        [pos[1][2],pos[1][2]])
    
    plt.show()