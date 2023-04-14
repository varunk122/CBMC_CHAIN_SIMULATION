from config import *
import os
import numpy as np

os.mkdir(project_name)

def write_positions(positions, step):
    os.mkdir(f"{project_name}/{step}")
    np.save(f'{project_name}/{step}/positions.npy', np.array(positions, dtype=object), allow_pickle=True)


