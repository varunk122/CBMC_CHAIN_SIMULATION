import numpy as np
from utility import pbc,add_bond_optimal,add_bond_with_torsion,add_bond_with_torsion_optimally

from config import *

def generate_random_atom():
    r = np.random.rand(3) * box_length
    return [r]

def generate_random_molecule(bond_length):
    r1 = np.zeros(3)
    r2 = np.random.rand(3)
    scale = bond_length / np.linalg.norm(r2)
    r2 = r2 * scale

    #translate

    r = np.random.rand(3)* box_length
    r1 += r
    r2 += r
    return [r1,r2]

def generate_propane_molecule(bond_length):
    mol_pos = generate_random_molecule(bond_length)
    r3 = add_bond_optimal(mol_pos)
    mol_pos.append(r3)
    return mol_pos

def generate_butane_molecule(bond_length):
    mol_pos = generate_propane_molecule(bond_length)
    r4 = add_bond_with_torsion(mol_pos)
    mol_pos.append(r4)
    return mol_pos

# def generate_random_propane_molecule(bond_length)
def check_overlap(mol_pos, positions):
    threshold  = bond_length*1.5
    for position in positions:
        for atom_pos in position:
            # print(np.linalg.norm(atom_pos - mol_pos[0]))
            for atom_mol in mol_pos:

                dx, dy, dz  = atom_pos - atom_mol
                dx, dy, dz = pbc(dx,box_length) , pbc(dy,box_length) , pbc(dz,box_length)
                r = np.sqrt(dx*dx+dy*dy+dz*dz)
            # dx, dy, dz  = atom_pos - mol_pos[1]
            # dx, dy, dz = pbc(dx,box_length) , pbc(dy,box_length) , pbc(dz,box_length)
            # r2 = np.sqrt(dx*dx+dy*dy+dz*dz)

                if r <threshold:
                    return True
    return False



def init_system(box_length, Npart):

    positions = []
    for i in range(Npart):
        max_iter = 1000
        it = 0
        while it < max_iter:
            if molecule_type == 'methane':
                mol_pos = generate_random_atom()
            elif molecule_type == 'ethane':
                mol_pos = generate_random_molecule(bond_length)
            elif molecule_type == 'propane':
                mol_pos = generate_propane_molecule(bond_length)
            elif molecule_type == 'butane':
                mol_pos = generate_butane_molecule(bond_length)
            else :
                print("Error: Unknown molcule_type in config file ")
                quit()
            if check_overlap(mol_pos, positions) == True:
                it += 1
            else :
                positions.append(mol_pos)
                break
        
        if it == max_iter:
            break

    assert len(positions) == Npart, f"Only {len(positions)} partcles are inserted instead of {Npart}"

    return positions


