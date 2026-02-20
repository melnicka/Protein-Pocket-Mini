import biotite.structure as struct
from biotite.structure import AtomArray, CellList
import numpy as np
from .config import Config


def find_pockets( protein_arr: AtomArray,
    ligand_list: list[AtomArray],
    cfg: Config
) -> list[AtomArray]:
    """Finds all of the protein's ligand binding pockets, based on the proximity to the ligand.

        Filters out solvent and stray ions, retaining polymers only.

        Args:
            protein_arr: AtomArray of the target protein.
            ligand_list: A list of AtomArrays, each representing a distinct ligand.
            cfg: Configuration object with an attribute:
                - radius (float): search radius in Angstroms [Å] used to identify binding pocket.

        Returns:
            A list of AtomArrays representing the filtered binding pockets. 
            If a search yields zero valid polymer atoms, that list index will contain None.
        """
    cell_list = CellList(protein_arr, cell_size=cfg.radius)
    pockets = []
    for ligand_arr in ligand_list:
        coords = ligand_arr.coord
        raw_indicies = cell_list.get_atoms(coords, radius=cfg.radius)
        indicies = np.unique(raw_indicies[raw_indicies != -1])

        raw_pocket = protein_arr[indicies]
        pocket_mask = (
            (struct.filter_polymer(raw_pocket)) |
            (struct.filter_polymer(raw_pocket, pol_type='nucleotide')) |
            (struct.filter_polymer(raw_pocket, pol_type='carbohydrate'))
        )
        pocket = raw_pocket[pocket_mask]
        pocket.bonds = struct.connect_via_residue_names(pocket)
        if len(pocket) == 0: 
            pocket = None 
        pockets.append(pocket)

    return pockets

def calculate_descriptors(pocket: AtomArray) -> dict:
    """Calculates descriptors of a binding pocket.

        Args:
            pocket: AtomArray of the pocket.

        Returns:
            A dictionary containing descriptors:
            - atom_count (int): Total number of atoms in the array.
            - gyration_radius (float): Spatial compactness of the pocket.
            - sasa (float): Total solvent accessible surface area.
            - Additional chemical composition metrics merged from get_res_info().
        """
    descriptors = {}
    descriptors['atom_count'] = len(pocket)
    descriptors['gyration_radius'] = float(struct.gyration_radius(pocket))

    sasa = struct.sasa(pocket)
    descriptors['sasa'] = float(np.sum(sasa))

    descriptors.update(get_res_info(pocket))

    return descriptors

def get_res_info(protein_arr: AtomArray):
    """Calculates chemical composition of the pocket by categorizing its atoms.

        Args:
            protein_arr: AtomArray of the pocket.

        Returns:
            A dictionary containing:
            - aromatic_count (int): Total number of individual atoms that belong 
              to aromatic residues (PHE, TYR, TRP, HIS).
            - hydrophobic_perc (float): Percentage of atoms belonging to nonpolar residues.
        """
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

    _, unique_indices = np.unique(protein_arr.res_id, return_index=True)
    unique_residues = protein_arr[unique_indices]

    for res in unique_residues:
        if res.res_name in aromatic:
            aromatic_count += 1
        if res.res_name in nonpolar:
            nonpolar_count += 1
        elif res.res_name in polar:
            polar_count += 1

    res_info['aromatic_count'] = aromatic_count
    res_info['hydrophobic_perc'] = 100 * nonpolar_count / (polar_count+nonpolar_count)

    return res_info
