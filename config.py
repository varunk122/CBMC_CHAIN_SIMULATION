#parameters
rho = 0.01      #density of the system
Npart = 108     #number of molecules
box_length = (Npart / rho) ** (1/3)  #box length of the system
bond_length = 1.5  #angstrom       #bond length of the system
write_interval = 100            #write interval

k = 30         #number of trial orientations
temp = 300 #kelvin     #temperature of the system
beta = 1               #value of beta / KbT

sigma = 3.75   #for LJ potential
eps = 0.46   # e/KbT #for LJ Potential
rcut = 3*bond_length #only for calculating pressure
theta_mean = 1.936 #radian
k_prop = 212 # KT/rad^2 

nsteps = 200000  #total number of simulation steps

project_name = "run1" #name of the project
molecule_type = "propane" #name of molecule to simulate
