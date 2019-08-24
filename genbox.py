#!/usr/bin/env python

import sys 
import os 
import gromacs
import numpy as np
solute = sys.argv[1]
outname = solute.split(".")[0]+"_b"+".gro"
gromacs.editconf(f=solute, box=list((input("Enter the box size: ").split())), o=outname)
# center=list(input("Enter the position for solute: ").split()),

