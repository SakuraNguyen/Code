#!/usr/bin/python 
import MDAnalysis as mda
import numpy as np
import sys
import math

u = mda.Universe(sys.argv[1], sys.argv[2])
#print(u.residues('Protein'))
#select the set of residues 

prot = u.select_atoms("protein")
mem = u.select_atoms("resname POPC")

#Final max: 54.015029511461194
#Final min: 6.19760425698905

bins = np.linspace(10.,50.,21)
#pos = np.zeros(21)
pos = [None]*len(bins)
currentd = np.ones(len(bins))*100
#print(len(id_f))
#print(currentd)
#maximum = 0.0
#minimum = u.dimensions[0]
dis = []

#CALCULATION
for i,ts in enumerate(u.trajectory[:430]):
    prot_com = prot.residues.center_of_mass()[:2]
    mem_com = mem.residues.center_of_mass()[:2]
    d = np.linalg.norm(prot_com-mem_com)

    for j,ele in enumerate(bins): # 50, 45, 40...
        e = abs(d-ele)
        if e < currentd[j]:
            pos[j] = np.copy(ts.positions)
            currentd[j] = e

# WRITING OF THE SELECTED FRAMES
for i,crd in enumerate(pos):
    u.atoms.positions=crd
    outname = 'conf'+str(i)+'.gro'
    u.atoms.write(outname)

    # JUST TO VERIFY ONCE
    if True:
        prot_com = prot.residues.center_of_mass()[:2]
        mem_com = mem.residues.center_of_mass()[:2]
        d = np.linalg.norm(prot_com-mem_com)
        print ('Distance in '+outname+' is',d)


























#   if d > maximum:
#       maximum = d
#   if d < minimum:
#       minimum = d
#print "Final max:", maximum
#print "Final min:", minimum
#print("Man je drza jako opice")
