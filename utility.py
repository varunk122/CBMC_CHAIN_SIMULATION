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

def lj_potential(atom_pos1, atom_pos2,box_length,i,j):
    dx, dy, dz  = atom_pos1 - atom_pos2
    dx, dy, dz = pbc(dx,box_length) , pbc(dy,box_length) , pbc(dz,box_length)
    r = dx*dx + dy*dy + dz * dz
    # print(r)
    #print(np.sqrt(r))
    fr6 = np.power(sigma[i][j]**2 / r,3)
    return 4*eps[i][j]*(fr6*(fr6-1))

def find_potential_between_alkanes(mol_pos1, mol_pos2,box_length):
    energy = 0
    i = 0
    for atom_pos_1 in mol_pos1:
        j = 0
        for atom_pos_2 in mol_pos2:
            energy += lj_potential(atom_pos_1, atom_pos_2,box_length,i,j)
            j = j+1
        i = i + 1
    
    return energy

def total_energy(positions,box_length):
    energy = 0
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            energy += find_potential_between_alkanes(positions[i], positions[j],box_length)

    return energy

def energy_of_a_chain(idx, positions, box_length):
    energy = 0
    for i in range(len(positions)):
        if i != idx:
            energy += find_potential_between_alkanes(positions[idx], positions[i], box_length)

    return energy

def energy_of_particle(idx1, idx2, positions,box_length):
    energy = 0
    for i in range(len(positions)):
        j = 0
        if i != idx1:
            for atom_pos in positions[i]:
                # r =  np.linalg.norm(positions[idx1][idx2] - atom_pos)
                energy += lj_potential(positions[idx1][idx2], atom_pos, box_length,idx2,j)
                j+=1
    # if idx2 == 2:
        # energy += lj_potential(positions[idx1][0],positions[idx1][2],box_length)
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

def calculate_pressure_for_monoatomic(atom_pos1, atom_pos2,i,j):
    Vol = box_length**3

    dx = atom_pos1[0] - atom_pos2[0] 
    dy = atom_pos1[1] - atom_pos2[1]
    dz = atom_pos1[2] - atom_pos2[2]
    xpbc = pbc(dx, box_length)
    ypbc = pbc(dy, box_length)
    zpbc = pbc(dz, box_length)
    r2 = xpbc*xpbc + ypbc*ypbc + zpbc*zpbc
    
    if r2 < rcut*rcut:                                              ### Calculating energy for particles inside cutoff distance
        fr2 = (sigma[i][j]*sigma[i][j])/r2
        fr6 = fr2*fr2*fr2
        
        return 48*eps[i][j]*((fr6*(fr6 - 0.5))) /(3*Vol)
    
    return 0

def calculate_pressure_for_chain_molecules(mol_pos1, mol_pos2):
    vir = 0
    i = 0
    for ap1 in mol_pos1:
        j = 0
        for ap2 in mol_pos2:
            vir += calculate_pressure_for_monoatomic(ap1, ap2,i,j)
            j += 1
        i += 1
    return vir

def virial_pressure_molecular(idx, positions):
    vir = 0
    for i in range(Npart):
        if i != idx:
            vir += calculate_pressure_for_chain_molecules(positions[idx], positions[i])
    return vir

def calculate_pressure(positions):
    frcut3 = (sigma[0][0]**3)/(rcut*rcut*rcut)
    Vir = 0
    Vol = box_length**3
    rho = Npart/Vol
    Ptail = (16*3.14/3)*rho*rho*eps[0][0]*(sigma[0][0]**3)*( (2*frcut3*frcut3*frcut3/3) - frcut3)      ### Tail correction to Pressure (Refer Frenkel and Smit for expression)
    #print(Ptail)
    for i in range(Npart):                                                ### Loop over all Virial interactions
        for j in range(i+1, Npart):
            Vir += calculate_pressure_for_chain_molecules(positions[i],positions[j])
    
    ### contribution of Ideal + Virial + Tail correction to pressure
    Pressure = (rho  / beta + Vir + Ptail)                
    return Pressure



def add_bond(mol_pos):
    r = mol_pos[-2] - mol_pos[-1]
    r = r / np.linalg.norm(r)
    flag = False
    while flag == False:
        r_new = generate_random_unit_vector()
        theta = np.arccos(np.dot(r,r_new))

        energy = 0.5*k_prop*(theta- theta_mean)**2
        # print(theta)
        # print(energy)
        
        if random.random() < np.exp(-beta*energy):
            r_new  = r_new *  bond_length + mol_pos[-1]
            flag = True
            return r_new

def convert_to_unit_vector(vector):
    return vector / np.linalg.norm(vector)

def add_bond_optimal(mol_pos):
    z = mol_pos[-2] - mol_pos[-1]
    y = np.array([0,1,0]) 
    if z[1] != 0:
        y = np.array([1,-(z[0]+z[2])/z[1],1])

    x = np.cross(y,z)

    #convert to unit vectors
    z = convert_to_unit_vector(z)
    y = convert_to_unit_vector(y)
    x = convert_to_unit_vector(x)

    #adding bond in spherical coordintate system for ease 

    r = bond_length
    theta = theta_mean
    phi = random.random() * (2*np.pi)

    #converting this back to cartesian coordintates

    return mol_pos[-1] + (r*np.sin(theta)*np.cos(phi))*x + (r*np.sin(theta)*np.sin(phi))*y + (r*np.cos(theta))*z

def add_bond_with_torsion(mol_pos):
    vjk = mol_pos[-1] - mol_pos[-2]
    vij = mol_pos[-2] - mol_pos[-3]
    
    c1 = np.cross(vij, vjk)
    while True:
        r4 = add_bond_optimal(mol_pos)
        vkl = r4 - mol_pos[-1]
        c2 = np.cross(vjk, vkl)
        phi = np.arccos(np.dot(c1,c2)/ (np.linalg.norm(c1)*np.linalg.norm(c2)))
        #parameters for butane in KJ/mol
        utors = 2.95*(1+np.cos(phi)) - 0.566*(1-np.cos(2*phi)) + 6.576*(1+np.cos(3*phi))

        if random.random() < np.exp(-beta * utors):
            # print(phi)
            return r4
        
def add_bond_with_torsion_optimally(mol_pos):
    vjk = mol_pos[-1] - mol_pos[-2]
    vij = mol_pos[-2] - mol_pos[-3]

    c1 = np.cross(vij, vjk)
    c2 = np.cross(c1,vjk)
    
    c1_ = convert_to_unit_vector(c1)
    c2_ = convert_to_unit_vector(c2)
    c3 = convert_to_unit_vector(vjk)

    phi = 0
    while True:
        r4 = add_bond_optimal(mol_pos)
        vkl = r4 - mol_pos[-1]
        c2 = np.cross(vjk, vkl)
        phi = np.arccos(np.dot(c1,c2)/ (np.linalg.norm(c1)*np.linalg.norm(c2)))
        #parameters for butane in KJ/mol
        utors = 2.95*(1+np.cos(phi)) - 0.566*(1-np.cos(2*phi)) + 6.576*(1+np.cos(3*phi))

        if random.random() < np.exp(-beta * utors):
            print(phi)
            v1 = mol_pos[-1]+ bond_length*np.sin(theta_mean)*(np.cos(phi)*c2_ + np.sin(phi)*c1_) + bond_length*np.cos(theta_mean)*c3
            v2 = mol_pos[-1]+ bond_length*np.sin(theta_mean)*(np.cos(phi)*c2_ - np.sin(phi)*c1_) + bond_length*np.cos(theta_mean)*c3
            print(v1,v2, r4)
            return r4 
        
    vkl = bond_length*np.sin(theta_mean)*(np.cos(phi)*c2 + np.sin(phi)*c1) + bond_length*np.cos(theta_mean)*c3

def read_positions_from_file(file_name):
    np_positions = np.load(file_name, allow_pickle=True)
    positions = []
    for mol_pos in np_positions:
        mol_pos_list = []
        for atom_pos in mol_pos:
            mol_pos_list.append(atom_pos)
        positions.append(mol_pos_list)
    
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
        f.write("%5i"%(3*N))
        f.write('\n')
        for i in range(N):
            f.write("%5i%.5s%5s%5i%8.3f%8.3f%8.3f"%(i+1,'ETH','C1',3*i+1, r[i][0][0], r[i][0][1], r[i][0][2]))
            f.write('\n')
            f.write("%5i%.5s%5s%5i%8.3f%8.3f%8.3f"%(i+1,'ETH','C2',3*i+2, r[i][1][0], r[i][1][1], r[i][1][2]))
            f.write('\n')
            f.write("%5i%.5s%5s%5i%8.3f%8.3f%8.3f"%(i+1,'ETH','C3',3*i+3, r[i][2][0], r[i][2][1], r[i][2][2]))
            f.write('\n')
        f.write("%10.5f%10.5f%10.5f"%(box_length,box_length,box_length))
        f.write('\n')
    