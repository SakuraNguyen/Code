#!/usr/bin/python


import MDAnalysis as mda
import numpy as np
import sys

u = mda.Universe(sys.argv[1], sys.argv[2])

prot = u.select_atoms("protein")
mem = u.select_atoms("resname POPC")

dis = []
outdist = open("z_distance.dat", "w")
print("CALCULATING.....")
#CALCULATION DISTANCE ALONG Z AXIS
for i,ts in enumerate(u.trajectory[:]):
#    print(i, ts.frame)
    prot_com = prot.residues.center_of_mass()[2]

    mem_com = mem.residues.center_of_mass()[2]
    d = float(prot_com) - float(mem_com)
    dis.append(d)
#    print(d)
    outdist.write('{:8d} {:10.2f}\n'.format(i, d))
outdist.close()
#    d = np.linalg.norm(prot_com-mem_com)
#print(dis)
#print(len(dis))
print(max(dis))
print(min(dis))
