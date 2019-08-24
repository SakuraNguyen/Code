#!/usr/bin/env python 
import os 
import sys

#path = os.getcwd()
#os.mkdir('combined_pull_file')
filename = [sys.argv[1], sys.argv[2]]

with open('tmp.xvg', 'w') as tmpfile:
    for line in open(sys.argv[2]):
        li=line.strip()
        if not (li.startswith("#") or li.startswith("@")):
            tmpfile.write(line)
newfilename = [sys.argv[1],"tmp.xvg"]
with open("out.xvg", "w") as outfile:
    for fname in newfilename:
        for lines in open(fname):
            outfile.write(lines)
outfile.close()





