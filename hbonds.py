#!/usr/bin/env python
import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree
from packing import make_whole
import itertools
import sys


def transformation(a, b, q, p):
    """
    Transform coordinates q when centering p (p=0)
    with periodic boundary conditions
    a... left border
    b... right border

    """
    return a + np.mod(q - p - a, b -a)

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
    if distOO < 3.5:
        H_atoms = [H11, H12, H21, H22]
        hb = 0
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
#                continue
#            elif angle > 120:
#                hb = hb + 1
    else:
        return 0

u = mda.Universe(sys.argv[1], sys.argv[2])
hbs_all_frame = np.empty((0,187,6), int)
for ts in u.trajectory[:]:
    box = u.dimensions[:3]
    u.atoms.positions = u.atoms.pack_into_box()
    Li_pos = u.select_atoms('name LI').positions
    #W = u.select_atoms('resname SOL')
    O_atom = u.select_atoms('name OW')
    O_pos = O_atom.positions
    kdtree = cKDTree(O_pos, boxsize=box)
    #print(len(kdtree.data))
    _, ind = kdtree.query(Li_pos, 4)
    #closest_neighbor = kdtree.query_ball_point(comw,2.5)
#    print(ind)
    hbs_one_frame = np.empty((0,6), int)
    for i in range(0, len(ind)):
    #    if i == 1:
    #        break
        id_o = ind[i]
        closest_w = O_atom.residues[id_o]
        closest_w.atoms.write("test.gro")
        pairs = list(itertools.combinations(closest_w, 2))
#        print(pairs[0][0])
        hbs_one_Li = np.zeros(len(pairs))
    #    print(hbs_one_Li)
        for j, w in enumerate(pairs):
           # if j == 5:
            #    break
#            print(j, w)
            hbs = hbond(w[0], w[1])
            hbs_one_Li[j] = hbs
    #        print(hbs_one_Li)
#        print(hbs_one_Li.reshape(1,6))
        hbs_one_frame = np.append(hbs_one_frame, hbs_one_Li.reshape(1,6), axis=0)
#    print(hbs_one_frame.shape)
    hbs_all_frame = np.append(hbs_all_frame, hbs_one_frame.reshape(1, 187, 6), axis=0)
print(hbs_all_frame.shape)
np.save("Hbond_Li4", hbs_all_frame)
