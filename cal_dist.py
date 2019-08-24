`python
import MDAnalysis as mda
import numpy as np
u = mda.Universe('topology.gro', 'trajectory.xtc')
# select first set of residues; in this example resids 21 - 27
residues1 = u.select_atoms('protein and resid 21:27').residues
# select second set of residues; in this example resids 38 - 42
residues2 = u.select_atoms('protein and resid 38:42').residues
# calculate distances between each CA in first group with each CA in second group
from MDAnalysis.lib.distance import distance_array
dists = distance_array(residues1.select_atoms('name CA').positions, residues2.select_atoms('name CA').positions)
# or, calculate distances between the COM of each residue in first group with COM of each residue in second group
residues1_com = np.array([res.center_of_mass() for res in residues1])
residues2_com = np.array([res.center_of_mass() for res in residues2])
dists = distance_array(residues1_com, residues2_com)
