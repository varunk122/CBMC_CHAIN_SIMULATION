import numpy as np
from utility import pbc
from config import *

def generate_random_molecule(init_pos,box_length):
    r = np.random.rand(3)* box_length #generating random point for center of mass
    scale=abs(r)
    #print(scale)
    #init_pos
    return init_pos-scale

def check_overlap(mol_pos, positions):
    threshold  = 0.3
    for position in positions:
        for atom_pos in position:
            # print(np.linalg.norm(atom_pos - mol_pos[0]))

            dx, dy, dz  = atom_pos - mol_pos[0]
            dx, dy, dz = pbc(dx,box_length) , pbc(dy,box_length) , pbc(dz,box_length)
            r1 = np.sqrt(dx*dx+dy*dy+dz*dz)
            dx, dy, dz  = atom_pos - mol_pos[1]
            dx, dy, dz = pbc(dx,box_length) , pbc(dy,box_length) , pbc(dz,box_length)
            r2 = np.sqrt(dx*dx+dy*dy+dz*dz)

            if r1 <threshold or r2 < threshold:
                return True
    return False

def init_system(init_pos,box_length, Npart):

    positions = []
    for i in range(Npart):
        max_iter = 1000
        it = 0
        while it < max_iter:
            r1,r2 = generate_random_molecule(init_pos,box_length)
            #print(r1)
            mol_pos = [r1,r2]
            if check_overlap(mol_pos, positions) == True:
                it += 1
            else :
                positions.append(mol_pos)
                #print(mol_pos)
                break
        
        if it == max_iter:
            break
    box_center=np.array([box_length/2,box_length/2,box_length/2])
    center=np.mean(np.mean(positions,axis=0),axis=0)
    return positions-center+box_center # translate to middle of box


