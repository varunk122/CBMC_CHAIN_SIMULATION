import math
import numpy as np
import matplotlib.pyplot as plt

#parameters

rho = 0.01
Npart = 100

box_length = (Npart / rho) ** (1/3)


def generate_random_molecule(bond_length, box_length):
    r1 = np.zeros(3)
    r2 = np.random.rand(3)
    scale = bond_length / np.linalg.norm(r2)
    r2 = r2 * scale

    #translate

    r = np.random.rand(3)* box_length
    r1 += r
    r2 += r
    return r1, r2

def check_overlap(mol_pos, positions,box_length):
    threshold  = 0.01
    for position in positions:
        for atom_pos in position:
            # print(np.linalg.norm(atom_pos - mol_pos[0]))
            if np.linalg.norm(atom_pos - mol_pos[0]) <threshold or np.linalg.norm(atom_pos - mol_pos[1]) < threshold:
                return True
    return False


def visualize(positions):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for pos in positions:
        ax.plot([pos[0][0],pos[1][0]],[pos[0][1],pos[1][1]],
        [pos[1][2],pos[1][2]])
    
    plt.show()

def init_system(box_length, Npart):

    positions = []
    for i in range(Npart):
        max_iter = 100
        it = 0
        while it < max_iter:
            r1,r2 = generate_random_molecule(1, box_length)
            mol_pos = [r1,r2]
            if check_overlap(mol_pos, positions, box_length) == True:
                it += 1
            else :
                positions.append(mol_pos)
                break
        
        if it == max_iter:
            break


    return positions


positions = init_system(box_length, Npart)
print(box_length)
visualize(positions)

