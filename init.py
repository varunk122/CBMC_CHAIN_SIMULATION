import numpy as np
from utility import pbc

from config import *

def generate_random_molecule(bond_length ):
    r1 = np.zeros(3)
    r2 = np.random.rand(3)
    scale = bond_length / np.linalg.norm(r2)
    r2 = r2 * scale

    #translate

    r = np.random.rand(3)* box_length
    r1 += r
    r2 += r
    return r1, r2

# def generate_random_propane_molecule(bond_length)
def check_overlap(mol_pos, positions):
    threshold  = 2
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



def init_system(box_length, Npart):

    positions = []
    for i in range(Npart):
        max_iter = 1000
        it = 0
        while it < max_iter:
            r1,r2 = generate_random_molecule(bond_length)
            mol_pos = [r1,r2]
            if check_overlap(mol_pos, positions) == True:
                it += 1
            else :
                positions.append(mol_pos)
                break
        
        if it == max_iter:
            break


    return positions


