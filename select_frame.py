#!/usr/bin/env python

"""The usage of this code
To run the code, four arguments should be provided in the order:
    1. distance between protein and membrane
    2. topology file and xtc file
    3. umbrella.mdp file
Be careful to change some parameters corresponding to the specific system,
for example the max, min distance, the window spacing, topology filename.
"""

import numpy as np
import os
import MDAnalysis as mda
import sys
import gromacs
import fileinput
import shutil

# Load the file containing the frame and the corresponding distance
data = np.loadtxt(sys.argv[1], usecols=1)
dist = []
index = []
for i in np.arange(30, 55, 2):
    
    substract = []
    for j, num in enumerate(data):
        val = abs(i-num)
        substract.append(val)
    #print(np.asarray(substract).shape)
    ID = substract.index(min(substract))
    index.append(ID)
    dist.append(data[ID])
print(index)
print(dist)

# Generate the frames for umbrella sampling
u = mda.Universe(sys.argv[2], sys.argv[3])
for ts in u.trajectory[np.array(index)[:]]:
    print ("Frame: {0:7d}, timestep: {1:12.4f}".format(ts.frame, ts.time))
    groname = "conf"+str(ts.frame)+".gro"
    u.atoms.write(groname)

# creating the dictionary for pairing the elements in index list and other in distance list as key:value
dic = dict(zip(index, dist))
for i, key in enumerate(list(dic.keys())):
    if i == 2:
        break
    print(i, key, dic[key])
    shutil.copy(sys.argv[5], "umbrella_cp.mdp")
    with fileinput.FileInput("umbrella_cp.mdp", inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace("CHANGE", str(dic[key]/10)), end='')
    gromacs.grompp(f="umbrella_cp.mdp", c="conf"+str(key)+".gro", p="topols.top", n="index.ndx", o="umbrella"+str(key)+".tpr")        
