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

def id2res(ind, group, molecule=True):
    """
    array_like ind : array of id of residues
    """
    mols = group[ind]
    if molecule:
        mols = mols.residues.atoms
    return mols

def neighbors2shells(ind,group,num_cut=4,molecule=True):
    first, second = np.array_split(ind,[4], axis=1)

    shell_pairs = []
    for i in range(0,len(ind)):
        sh1 = id2res(first[i],group, molecule)
        sh2 = id2res(second[i],group, molecule)
        shell_pairs.append([sh1,sh2])
    return shell_pairs



DEBUG = False

if __name__ == '__main__':

    u = mda.Universe(sys.argv[1], sys.argv[2])

    all_c1 = 0
    all_num_hb = 0

    for ts in u.trajectory[:1]:

        box = u.dimensions[:3]

        u.atoms.positions = u.atoms.pack_into_box()
        Li_pos = u.select_atoms('name LI').positions
        O_atoms= u.select_atoms('name OW')
        O_pos = O_atoms.positions

        kdtree = cKDTree(O_pos, boxsize=box)
        _, ind = kdtree.query(Li_pos, 16)
        shell_pairs = neighbors2shells(ind,O_atoms)

        c1 = 0
        num_hb = 0
        for pairs in shell_pairs:

            sh1 = pairs[0]
            sh2 = pairs[1]

            if DEBUG:
                print ("SHELL")
            for elem_sh1 in sh1.residues:
                for elem_sh2 in sh2.residues:
                    num_hb += hbond(elem_sh1,elem_sh2)
        ave_one_frame = num_hb/187
        print(ave_one_frame)
#                    if DEBUG:
#                        print (elem_sh1,elem_sh2,hbond(elem_sh1,elem_sh2))
    #                 c1 =+ 1
    #     all_c1 += c1
    #     all_num_hb += num_hb
    #
    #     print(float(num_hb)/c1)
    #
    # print('OVERALL AVERAGE HBONDS AROUND LITHIUMS IN THE SYSTEM:', float(all_num_hb)/all_c1)
