import biotite.structure as struct
from biotite.structure import AtomArray, CellList
import numpy as np
from .config import Config

def find_pockets(
    protein_arr: AtomArray,
    ligand_list: list[AtomArray],
    cfg: Config
) -> list[AtomArray]:
    cell_list = CellList(protein_arr, cell_size=cfg.radius)
    pockets = []
    for ligand_arr in ligand_list:
        coords = ligand_arr.coord
        raw_indicies = cell_list.get_atoms(coords, radius=cfg.radius)
        indicies = raw_indicies[raw_indicies != -1]

        raw_pocket = protein_arr[indicies]
        pocket_mask = (
            (struct.filter_polymer(raw_pocket)) |
            (struct.filter_polymer(raw_pocket, pol_type='nucleotide')) |
            (struct.filter_polymer(raw_pocket, pol_type='carbohydrate'))
        )
        pocket = raw_pocket[pocket_mask]
        pocket.bonds = struct.connect_via_residue_names(pocket)
        pockets.append(pocket)

    return pockets

def calculate_descriptors(pocket: AtomArray) -> dict:
    descriptors = {}
    descriptors['atom_count'] = len(pocket)
    descriptors['gyration_radius'] = float(struct.gyration_radius(pocket))

    sasa = struct.sasa(pocket)
    descriptors['sasa'] = float(np.sum(sasa))

    descriptors.update(get_res_info(pocket))

    return descriptors

def get_res_info(protein_arr: AtomArray):
    res_info = {}
    nonpolar = [
        "GLY", "ALA", "VAL", "LEU", "ILE",
        "MET", "PHE", "TRP", "PRO"
    ]
    polar = [
        "SER", "THR", "ASN", "GLN", "CYS", "TYR",
        "ASP", "GLU", "LYS", "ARG", "HIS"
    ]
    aromatic = ["PHE", "TYR", "TRP", "HIS"]
    nonpolar_count, polar_count, aromatic_count = 0, 0, 0

    for atom in protein_arr:
        if atom.res_name in aromatic:
            aromatic_count += 1
        if atom.res_name in nonpolar:
            nonpolar_count += 1
        elif atom.res_name in polar:
            polar_count += 1

    res_info['aromatic_count'] = aromatic_count
    res_info['hydrophobic_perc'] = 100 * nonpolar_count / (polar_count+nonpolar_count)

    return res_info
