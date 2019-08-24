#!/usr/bin/env python 
import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree
from packing import make_whole


def distanceOO(w1, w2):
    """
    Calculate the distance between two water molecules through their O atoms
    """
    O1 = 
def hbond_angle(w1,w2):


    # minimum image!
    make_whole(w1)
    make_whole(w2)

    o1  = w1.atoms[0].position
    h11 = w1.atoms[1].position
    h12 = w1.atoms[2].position
    o2  = w2.atoms[0].position
    h21 = w2.atoms[1].position
    h22 = w2.atoms[2].position
    hs = [h11, h12, h21, h22]
    for h in hs:
        v1 = o1-h
        v2 = o2-h
        lv1 = np.linalg.norm(v1)
        lv2 = np.linalg.norm(v2)
        angle = np.arccos(np.clip(np.dot(v1, v2)/(lv1*lv2), -1.0, 1.0))
        angle *= 180./np.pi
        if angle > angle_oho:
            return 1
    return 0

ncut = 2.5
doo = 3.5
angle_oho = 120.

u = mda.Universe('03-water-LIsCLs-pro-mol.gro')
box = u.dimensions[:3]
Li = u.select_atoms('name LI')
Wa = u.select_atoms('resname SOL')
hbonds = 0
waters = 0

# THIS NEEDS TO LOOP OVER THE TRAJ

u.atoms.positions = u.atoms.pack_into_box()
comWa = Wa.select_atoms('name OW').positions
comLi = Li.positions
kdtree = cKDTree(comWa, boxsize=box)

Li_nlist = kdtree.query_ball_point(comLi,ncut)

for i in range(0,len(Li_nlist)):  # looping over the balls of water around Li atoms
    nbs = Li_nlist[i]
    w = Wa.atoms.residues[nbs]
    comw = w.atoms.select_atoms('name OW').positions
    W_nlist = kdtree.query_ball_point(comw,doo)

    for j in range(0,len(W_nlist)):  # looping over the waters in a ball
        waters += 1
        wnbs = W_nlist[j]
        idx = w[j].resid
        close_waters = Wa.atoms.residues[wnbs]

        for k in close_waters:
            hbonds += hbond_angle(w,k)
            
print (float(hbonds)/waters)
