import init
import utility
import visualize
import random
import numpy as np
import matplotlib.pyplot as plt

#import parameters
from config import *

debug = False
accepted_steps = 0


def CBMC_step(positions, Npart):
    #select a random_chain
    global accepted_steps
    idx = random.randint(0,Npart-1)
    #find rosenbluth factor for previous configuration
    Uo = utility.energy_of_a_chain(idx, positions, box_length)
    Po = utility.virial_pressure_molecular(idx, positions)
    prev_position = np.copy(positions[idx])
    Wo1 = k * np.exp(-beta * utility.energy_of_particle(idx,0,positions,box_length))
    Wo2 = np.exp(-beta * utility.energy_of_particle(idx,1,positions,box_length))
    for i in range(k-1):
        r = utility.generate_random_unit_vector() * bond_length
        positions[idx][1] = positions[idx][0] + r
        Wo2 += np.exp(-beta * utility.energy_of_particle(idx,1,positions,box_length))
    
    positions[idx][1] = prev_position[1]
    Wo3 = np.exp(-beta * utility.energy_of_particle(idx,2,positions,box_length))
    mol_pos = [prev_position[0], prev_position[1]]
    for i in range(k-1):
        r3 = utility.add_bond_optimal(mol_pos)
        positions[idx][2] = r3
        Wo3 += np.exp(-beta*utility.energy_of_particle(idx,2,positions,box_length))

    positions[idx][2] = prev_position[2]

    Wo4 = np.exp(-beta* utility.energy_of_particle(idx,3,positions,box_length))
    mol_pos = [prev_position[0], prev_position[1], prev_position[2]]
    for i in range(k-1):
        r4 = utility.add_bond_with_torsion(mol_pos)
        positions[idx][3] = r4
        Wo4 += np.exp(-beta * utility.energy_of_particle(idx,3,positions,box_length))

    positions[idx][3] = prev_position[3]

    Wo = Wo1 * Wo2 * Wo3 * Wo4

    #find rosenbluth factor for new configuration
    #choose a random postion for first atom
    positions[idx][0] = np.random.rand(3)* box_length
    Wn1 =  k * np.exp(-beta * utility.energy_of_particle(idx,0,positions,box_length))
    second_atom_pos = []
    Wn2 = []
    for i in range(k):
        r = utility.generate_random_unit_vector()* bond_length
        positions[idx][1] = positions[idx][0] + r
        second_atom_pos.append(positions[idx][1].copy())
        Wn2.append(np.exp(-beta * utility.energy_of_particle(idx,1,positions,box_length)))

    Wn2_sum = sum(Wn2)

    #select a configuration i with probaility Wn2[i]/Wn2_sum
    cum_Wn2 = Wn2[0]
    r_Wn2_sum = random.random()*Wn2_sum
    i = 0
    while cum_Wn2 < r_Wn2_sum:
        i +=1
        cum_Wn2 += Wn2[i]
    # Wn = Wn1 * Wn2_sum
    #replace position of second atom with selected second atom configuration 
    positions[idx][1] = second_atom_pos[i].copy()

    third_atom_pos = []
    Wn3 = []
    mol_pos = [positions[idx][0],positions[idx][1]]
    for i in range(k):
        r3 = utility.add_bond_optimal(mol_pos)
        positions[idx][2] = r3.copy()
        third_atom_pos.append(r3.copy())
        Wn3.append(np.exp(-beta*utility.energy_of_particle(idx,2,positions,box_length)))

    Wn3_sum = sum(Wn3)
    #select a configuration i with probability Wn3[i]/Wn3_sum
    cum_Wn3 = Wn3[0]
    r_Wn3_sum = random.random()* Wn3_sum
    i = 0
    while cum_Wn3 < r_Wn3_sum:
        i += 1
        cum_Wn3 += Wn3[i]
    
    positions[idx][2] = third_atom_pos[i].copy()

    fourth_atom_pos = []
    Wn4 = []
    mol_pos = [positions[idx][0], positions[idx][1], positions[idx][2]]
    for i in range(k):
        r4 = utility.add_bond_with_torsion(mol_pos)
        positions[idx][3] = r4.copy()
        fourth_atom_pos.append(r4.copy())
        Wn4.append(np.exp(-beta*utility.energy_of_particle(idx,3,positions,box_length)))
    
    Wn4_sum = sum(Wn4)
    #select a configuration i with probability Wn4[i]/Wn4_sum
    cum_Wn4 = Wn4[0]
    r_Wn4_sum = random.random() * Wn4_sum
    i = 0
    while cum_Wn4 < r_Wn4_sum:
        i+=1
        cum_Wn4 += Wn4[i]
    
    positions[idx][3] = fourth_atom_pos[i].copy()

    Wn = Wn1 * Wn2_sum * Wn3_sum * Wn4_sum

    if (Wn < Wo and random.random() > Wn/Wo) or (Wo == 0 and Wn==0) :
        #not accept
        positions[idx] = prev_position   
        print(f"Old Rosenbluth factor {Wo} New Rosenbluth factor {Wn}, Move not accpeted !! ")
        return positions, accepted_steps,0, 0
    else :
        accepted_steps += 1
        print(f"Old Rosenbluth factor {Wo} New Rosenbluth factor {Wn}, Move accpeted !! ")
        Un = utility.energy_of_a_chain(idx, positions, box_length)
        Pn = utility.virial_pressure_molecular(idx, positions)
        return positions, accepted_steps, Un - Uo, Pn - Po
