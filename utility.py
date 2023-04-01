import numpy as np
import random
import math

#parameters 

sigma = 1
eps = 0.466 # e/KbT


def pbc(x, L):
    # pbc checks for periodic image and returns nearest image
    # x is magnitude of distance between two atoms
    # L is box_length

    x = x - math.floor(x/L)* L

    if x >= (0.5*L):
        return L - x 
    
    return x

def lj_potential(r,box_length):
    r = pbc(r,box_length)
    # print(r)
    fr6 = sigma / np.power(r,6)
    return 4*eps*(fr6*(fr6-1))

def find_potential_between_alkanes(mol_pos1, mol_pos2,box_length):
    energy = 0
    for atom_pos_1 in mol_pos1:
        for atom_pos_2 in mol_pos2:
            r = np.linalg.norm(atom_pos_1 - atom_pos_2)
            energy += lj_potential(r,box_length)
    
    return energy

def total_energy(positions,box_length):
    energy = 0
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            energy += find_potential_between_alkanes(positions[i], positions[j],box_length)

    return energy

def energy_of_particle(idx1, idx2, positions,box_length):
    energy = 0
    for i in range(len(positions)):
        if i != idx1:
            for atom_pos in positions[i]:
                r =  np.linalg.norm(positions[idx1][idx2] - atom_pos)
                energy += lj_potential(r,box_length)
    
    # for i in range(idx2):
    #     r = np.linalg.norm(positions[idx1][i] - positions[idx1][idx2])
    #     energy += lj_potential(r,box_length)

    # print(energy)
    return energy

def generate_random_unit_vector():
    r = np.random.rand(3)
    scale = 1 / np.linalg.norm(r)
    r = r * scale

    return r

def get_minimum_distance(positions, idx1, idx2, box_length):
    min_len = 1000
    for i in range(len(positions)):
        if idx1 !=i:
            for pos in positions[i]:
                min_len = min(min_len, pbc(np.linalg.norm(pos - positions[idx1][idx2]),box_length))
    
    return min_len

def check(positions, bond_length):
    for pos in positions:
        if np.linalg.norm(pos[0] - pos[1]) > bond_length + 0.01:
            return False
    
    return True


# def add_bond(mol_pos, theta_mean, k):
#     r = abs(mol_pos[-2] - mol_pos[-1])
#     r = r / np.linalg.norm(r)
#     r_new = generate_random_unit_vector()
#     theta = np.arccos(np.dot(r,r_new))

#     energy = 0.5*k*(theta- theta_mean)^2
    
#     if random.random() < np.exp(-beta*energy):
#         r_new  = r_new + mol_pos[-1]
#         return r_new
    
#     return add_bond(mol_pos,theta_mean,k)
    