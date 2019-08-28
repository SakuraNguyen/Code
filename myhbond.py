#!/usr/bin/env python
import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree
from packing import make_whole
import itertools
import sys


def hbond(w1, w2):
    # minimum image!
    make_whole(w1)
    make_whole(w2)
    """
    Calculate the distance between two water molecules through their O atoms
    """
    O1 = w1.atoms[0].position
    H11 = w1.atoms[1].position
    H12 = w1.atoms[2].position
    O2 = w2.atoms[0].position
    H21 = w2.atoms[1].position
    H22 = w2.atoms[2].position
    distOO = np.linalg.norm(O1-O2)# bond length cutoff is of 3.0 Angstrom

    hb = 0
    if distOO < 3.5:
        H_atoms = [H11, H12, H21, H22]
        for H in H_atoms:
            v1 = O1-H
            v2 = O2-H
            lv1 = np.linalg.norm(v1)
            lv2 = np.linalg.norm(v2)
            angle = np.arccos(np.clip(np.dot(v1, v2)/(lv1*lv2), -1.0, 1.0))
            angle *= 180./np.pi
            if angle > 120:
                hb +=  1
    return hb

def pair_combination(*arrays):
    """
    This function is to generate pairs which are not duplicate from two arrays
    There are multiple parameters can input
    """
    for t in itertools.combinations(arrays, 2):
        for pair in itertools.product(*t):
            yield pair

if __name__ == '__main__':
    u = mda.Universe(sys.argv[1], sys.argv[2])
    for ts in u.trajectory[:1]:
        box = u.dimensions[:3]
        u.atoms.positions = u.atoms.pack_into_box()
        Li_pos = u.select_atoms('name LI').positions
        O_atom = u.select_atoms('name OW')
        O_pos = O_atom.positions
        kdtree = cKDTree(O_pos, boxsize=box)
        _, ind = kdtree.query(Li_pos, 16)
    #### Split each 16-element arry into two subarrays corresponding to the first
    #### and the second shell
        first, second = np.array_split(ind,[4], axis=1)
        print(len(second))
        first_res = []

        for i in range(0, len(first)): # loop over the number of Li atoms
    
            if i == 2:
                break
            id_o_first = first[i]
            closest_w_first = O_atom.residues[id_o_first]
            first_res.append(closest_w_first)
    #    print(first_res)
        second_res = []
        for j in range(0, len(second)): # loop over the number of Li atoms
            if j == 2:
                break
            id_o_second = second[j]
            closest_w_second = O_atom.residues[id_o_second]
            second_res.append(closest_w_second)
            pairs = []
        for n, p in enumerate(zip(first_res, second_res)):
            if n == 2:
                break
            pairs.append([x for x in pair_combination(p[0],p[1])])
        print(len(pairs[1]))
        for k, pa in enumerate(pairs):
            if k == 2:
                break
            for h, pai in enumerate(pa):
                # if h == :
                #     break
                hbs1 = hbond(pai[0], pai[1])
                print(hbs1)
