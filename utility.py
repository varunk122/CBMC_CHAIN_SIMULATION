import numpy as np
import random
import math

#parameters
from config import *


def pbc(dx, L):
    # pbc checks for periodic image and returns nearest image
    # x is magnitude of distance between two atoms
    # L is box_length
    if dx > L/2:
        dx -= L
    if dx < -L/2:
        dx += L
    
    return dx

def lj_potential(atom_pos1, atom_pos2,box_length):
    dx, dy, dz  = atom_pos1 - atom_pos2
    dx, dy, dz = pbc(dx,box_length) , pbc(dy,box_length) , pbc(dz,box_length)
    r = dx*dx + dy*dy + dz * dz
    # print(r)
    #print(np.sqrt(r))
    fr6 = np.power(sigma**2 / r,3)
    return 4*eps*(fr6*(fr6-1))

def find_potential_between_alkanes(mol_pos1, mol_pos2,box_length):
    energy = 0
    for atom_pos_1 in mol_pos1:
        for atom_pos_2 in mol_pos2:
            energy += lj_potential(atom_pos_1, atom_pos_2,box_length)
    
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
                # r =  np.linalg.norm(positions[idx1][idx2] - atom_pos)
                energy += lj_potential(positions[idx1][idx2], atom_pos, box_length)

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


def add_bond(mol_pos, theta_mean, k, beta = 1):
    r = abs(mol_pos[-2] - mol_pos[-1])
    r = r / np.linalg.norm(r)
    r_new = generate_random_unit_vector()
    theta = np.arccos(np.dot(r,r_new))

    energy = 0.5*k*(theta- theta_mean)^2
    
    if random.random() < np.exp(-beta*energy):
        r_new  = r_new + mol_pos[-1]
        return r_new
    
    return add_bond(mol_pos,theta_mean,k)

def read_positions_from_file(file_name):
    np_positions = np.load(file_name, allow_pickle=True)
    positions = []
    for mol_pos in np_positions:
        positions.append([mol_pos[0],mol_pos[1]])
    
    return positions

# def print_bond_length_of_mol(positions):

def read_positions_from_gromacs_file(file_name,Npart):
    with open(file_name,'r') as f:
        lines = [line for line in f]
        positions = []
        for i in range(Npart):
            l1 = lines[2*i+2].split()
            l2 = lines[2*i+3].split()
            mol_pos_1 = np.array([l1[3],l1[4],l1[5]]).astype(float)
            mol_pos_2 = np.array([l2[3],l2[4],l2[5]]).astype(float)
            positions.append([mol_pos_1,mol_pos_2])

    return positions

def output_energy(s, PE,file):
    with open(file, 'a+') as f:
        f.write("%9d %8.3f" %(s, PE))
        f.write('\n')

def output_xyz(N,box_length,r,p,file):
    """xyz output"""
    with open(file, 'a+') as f:
        f.write('step=%5d\n'%(p))
        f.write("%5i"%(2*N))
        f.write('\n')
        for i in range(N):
            f.write("%5i%.5s%5s%5i%8.3f%8.3f%8.3f"%(i+1,'ETH','C1',2*i+1, r[i][0][0], r[i][0][1], r[i][0][2]))
            f.write('\n')
            f.write("%5i%.5s%5s%5i%8.3f%8.3f%8.3f"%(i+1,'ETH','C2',2*i+2, r[i][1][0], r[i][1][1], r[i][1][2]))
            
            f.write('\n')
        f.write("%10.5f%10.5f%10.5f"%(box_length,box_length,box_length))
        f.write('\n')
    