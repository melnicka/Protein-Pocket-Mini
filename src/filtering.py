import biotite.structure as struct 
from .utils import get_atom_array

def find_ligand(atom_array: struct.AtomArray):
    ligands = atom_array[~struct.filter_polymer(atom_array)]
    ligands = ligands[~struct.filter_solvent(ligands)]

    ion_mask = struct.filter_monoatomic_ions(ligands)
    heavy_mask = struct.filter_heavy(ligands)
    rest_mask = (~ion_mask) & (~heavy_mask)
    ions = ligands[ion_mask]
    heavy = ligands[heavy_mask]
    rest = ligands[rest_mask] 

    return (ligands, ions, heavy, rest) 
