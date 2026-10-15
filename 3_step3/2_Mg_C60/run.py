
import os, time, copy, sys
import numpy as np
import scipy.sparse as sp

from libra_py import units
from libra_py import data_conv
import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.citools.ci as ci
from liblibra_core import *



labels, q = CP2K_methods.read_trajectory_xyz_file("C60_Mg-pos-1.xyz", 0)
print(labels)
print(q)

"""
From step2 setup, we have:

# Lowest and highest orbital, Here HOMO is 90
params['lowest_orbital'] = 90
params['highest_orbital'] = 160

so the lowest orbital is 90
"""

params = {"atom_labels":labels, "timestep":0, 
          "dt":1.0*units.fs2au,
          "nelec_act_space":None,
          "ci_threshold":0.001,
          "lowest_orbital": 90,
          "nstates": 11,
          "is_first_time":{0:True},
          "act_state":{0:1},
         }
print(params)

# Emulates 1 trajectory
full_id = Py2Cpp_int([0, 0])

res = "Mg_C60"
# Create working directory, if doesn't exist
if not os.path.exists(res):
    os.mkdir(res)

# Do the first 300 steps 
for i in range(300):
    #print(F"======== Iteration {i} ==============")
    labels, q = CP2K_methods.read_trajectory_xyz_file("C60_Mg-pos-1.xyz", i)
    params["timestep"] = i
    params["logfile_name"] = F"all_logfiles/step_{i}.log"
    params["time_overlap_filename"] = F"res/St_ks_{i}.npz"
    params["overlap_filename"] = F"res/S_ks_{i}.npz"
    
    obj = CP2K_methods.cp2k_compute_adi(q, params, full_id)        
    obj.ham_adi.show_matrix(F"{res}/ham_adi_{i}.txt")
    obj.hvib_adi.show_matrix(F"{res}/hvib_adi_{i}.txt")
    obj.time_overlap_adi.real().show_matrix(F"{res}/st_adi_{i}.txt")
    obj.overlap_adi.real().show_matrix(F"{res}/s_adi_{i}.txt")


