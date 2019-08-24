#!/usr/bin/env python

import gromacs
import sys, os

grofile = sys.argv[1]

gromacs.editconf(f=grofile, translate=list(input("Enter the translate vector: ").split()), o="out.gro")
