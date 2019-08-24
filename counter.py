#!/usr/bin/env python 
import sys
import os
from collections import Counter

list_0=[]
with open(sys.argv[1], "r") as f:
    for line in f.readlines():
        if 'POPS' in line:
            list_0.append(line.split()[0])
    print(list_0)
    counter=Counter(list_0)
    print(counter)
f.close()            


