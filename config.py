#parameters
rho = 0.015861948     #density of the system
Npart = 108     #number of molecules
box_length = (Npart / rho) ** (1/3)  #box length of the system
bond_length = 1.5  #angstrom       #bond length of the system
write_interval = 100            #write interval

k = 30         #number of trial orientations
temp = 100 #kelvin     #temperature of the system
beta = 1               #value of beta / KbT

sigma = 3.75   #for LJ potential
eps = 1.38   # e/KbT #for LJ Potential
rcut = 3*bond_length #only for calculating pressure
theta_mean = 1.936 #radian
k_prop = 212 # KT/rad^2 

nsteps = 1000000  #total number of simulation steps

project_name = "mc5" #name of the project
molecule_type = "methane" #name of molecule to simulate
