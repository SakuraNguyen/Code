#!/usr/bin/env python

import os 
import sys
import MDAnalysis as mda
import gromacs
import numpy as np

# LOAD THE GRO FILE OF MEMBRANE
u1 = mda.Universe(sys.argv[1])
tr_vector = [0, 0, u1.dimensions[2]/10]
#tr_input = (*tr_vector, sep = ' ')
#tr_vector= np.asarray([0, 0, u1.dimensions[2]/10]) 
#print(tr_vector)

###SET UP THE PROTEIN BOX 
print("This stage is to setup protein box to combine to membrane box!!!")
#boxsize = input("enter the box size: ").split()
z_box = float(input("enter the box size along z: "))
x_box = u1.dimensions[0]/10
y_box = u1.dimensions[1]/10
boxsize =[x_box, y_box, z_box]
#print(boxsize)
#print(type(z_box))
center = [x_box/2, y_box/2, float(input("enter the coordinates you want to place protein along z: "))]
filename = input("provide the file name of residue: ")
outname = filename.split(".")[0]
outname_b = outname+"_b"+".gro"
gromacs.editconf(f=filename, o=outname_b, box=boxsize, center=center)
outname_bw = outname+"_bw"+".gro"
gromacs.solvate(cp=outname_b, o="out.gro")
gromacs.editconf(f="out.gro", translate= tr_vector, o=outname_bw)
#
##
#### COMBINE MEMBRANE GRO FILE AND PEPTIDE GRO FILE
u2 = mda.Universe(outname_bw)
u12 = mda.Merge(u2.atoms, u1.atoms)

u12.atoms.write(str(sys.argv[1]).split(".")[0]+"_"+outname+".gro")
