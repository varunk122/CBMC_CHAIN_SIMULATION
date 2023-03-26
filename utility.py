import numpy as np
import math

#parameters 

sigma = 1
eps = 1

def pbc(x, L):
    # pbc checks for periodic image and returns nearest image
    if x >= (0.5*L):
        dist = x - L
    elif x < (-0.5*L):
        dist = x + L
    else:
        dist = x
    return dist

def lj_potential(r):
    r = pbc(r)
    fr6 = sigma / np.power(r,6)
    return 4*eps*(fr6*(fr6-1))

def total_energy(positions):
    energy = 0
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            r =  abs(positions[j] - positions[i])
            energy += lj_potential(r)

    return energy

def energy_of_particle(idx, positions):
    energy = 0
    for i in range(len(positions)):
        if i != idx:
            r =  abs(positions[idx] - positions[i])
            energy += lj_potential(r)

    return energy