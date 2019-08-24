import numpy as np
import MDAnalysis as mda

def pack_around(atom_group, center):
    """
    Translate atoms to their periodic image closest to a given point.

    The function assumes that the center is in the main periodic image.
    """
    assert(isinstance(atom_group,mda.core.groups.AtomGroup))
    box = atom_group.universe.dimensions
    atom_group.pack_into_box()
    positions = atom_group.positions.copy()
    sub = positions - center
    mask = np.where(np.sqrt(sub**2) > box[:3] / 2)
    positions[mask] -= (atom_group.universe.dimensions[mask[1]]  
                            * np.sign(sub[mask]))
    atom_group.positions = positions

def make_whole(residue):  
    """  
    Translate atoms so they are clustered around their COM.   
                           
    Only works for small molecules (relative to box)  
    """
    pack_around(residue.atoms, residue.atoms.center_of_mass())
