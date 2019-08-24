#!/usr/bin/env python

import MDAnalysis as mda
#from MDAnalysis.lib.distance import distance_array
import numpy as np
import os 
import sys


u = mda.Universe(sys.argv[1])
prot = u.select_atoms('protein')
#memb = u.select_atoms('resname POPC')
memb = u.select_atoms('resname POPE') + u.select_atoms('resname POPS')
print(prot)
print(memb)
prot_com = prot.residues.center_of_mass()
memb_com = memb.residues.center_of_mass()
#d = distance_array(prot_com,memb_com)
d = np.linalg.norm(prot_com-memb_com)
print(d)
