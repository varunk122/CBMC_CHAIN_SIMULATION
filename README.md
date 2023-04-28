# CBMC_CHAIN_SIMULATION

# Command to run the code using shell script

- Add these two lines in your bashrc 

PATH_CHAINMC=<absoulte_path_to_clone>

PATH=$PATH:$PATH_CHAINMC

Example - 
PATH_CHAINMC="/home/varun/Desktop/CBMC_CHAIN_SIMULATION"

PATH=$PATH:$PATH_CHAINMC

- And, also modify the PATH_CHAINMC variable in line 46 of chainmc file in the CBMC_CHAIN_SIMULATION folder as done in bashrc
- To Run the simulation use the following command

`chainmc.sh -c relative/path/to/config_file -l relative/path/to/log_file -t type`

Example Command - 

-  For MC

`chainmc - c ./config.py - l .  -t mc`

- For CBMC 

`chainmc - c ./config.py - l .  -t mc`

# Run the code without shell script 

- Go inside your CBMC_CHAIN_SIMULATION folder
- Modify config.py accordingly (currently_its_for_butane)
- Run `python3 cbmc_for_chains.py` for cbmc simulation and `python3 mc_for_chains.py` for mc respectively.

# LJ Potential Paramerters for different molecules
# Ethane
sigma = [[0.375,0.375],[0.375,0.375]]   #for LJ potential
eps =  [[0.814,0.814],[0.814,0.814]]  # KJ/mol #for LJ Potential
# Propane

sigma = [[0.375,0.385,0.375],[0.385,0.395,0.385],[0.375,0.385,0.375]]   #for LJ potential
eps =  [[0.814,0.558,0.814],[0.558,0.382,0.558],[0.814,0.558,0.814]]  # KJ/mol #for LJ Potential

# Butane

sigma = [[0.375,0.385,0.385,0.375],[0.385,0.395,0.395,0.385],[0.385,0.395,0.395,0.385],[0.375,0.385,0.385,0.375]]   #for LJ potential
eps =  [[0.814,0.558,0.558,0.814],[0.558,0.382,0.382,0.558],[0.558,0.382,0.382,0.558],[0.814,0.558,0.558,0.814]]  # KJ/mol #for LJ Potential

# Sample config file for Butane

```
#parameters
rho =   0.602  #density of the system
Npart = 108     #number of molecules
box_length =  (Npart / rho) ** (1/3)  #box length of the system
bond_length = 0.15  #nanometer       #bond length of the system
write_interval = 100            #write interval
k = 10         #number of trial orientations
temp = 180     #kelvin     #temperature of the system
beta = 1/(8.314*0.001 * temp)              #value of beta / KbT
sigma = [[0.375,0.385,0.385,0.375],[0.385,0.395,0.395,0.385],[0.385,0.395,0.395,0.385],[0.375,0.385,0.385,0.375]]   #for LJ potential
eps =  [[0.814,0.558,0.558,0.814],[0.558,0.382,0.382,0.558],[0.558,0.382,0.382,0.558],[0.814,0.558,0.558,0.814]]  # e/KbT #for LJ Potential
rcut = 3*0.375 #only for calculating pressure
theta_mean = 1.91 #radian
k_prop = 212 # KT/rad^2 

nsteps = 100000  #total number of simulation steps

project_name = "cbmc_for_butane" #name of the project
molecule_type = "butane" #name of molecule to simulate
```

