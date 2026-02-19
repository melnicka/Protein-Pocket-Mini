from biotite.structure.io import pdbx
from .data_fetching import get_protein_data
from .pocket import find_pockets, calculate_descriptors
from .config import Config

def display_full_info(pdb_id: str, cfg: Config):
    pdb_id = pdb_id.upper()
    protein_arr, ligands, metadata = get_protein_data(pdb_id, cfg)
    num_pockets = len(ligands)
    pockets = find_pockets(protein_arr, ligands, cfg)
    print(f"\nFOUND {num_pockets} POCKETS IN PROTEIN {pdb_id}")
    
    for i in range(0, num_pockets):
        print(f"\n-------------- POCKET {i+1} -------------------")
        pocket = pockets[i]
        descriptors = calculate_descriptors(pocket)
        lig_data = metadata[i]
        display_pocket_info(descriptors, lig_data)

    print("")

def display_pocket_info(pocket_descriptors: dict, metadata: dict):
    print(f"""------ LIGAND INFO ------
{metadata['name']}
Formula: {metadata['formula']}
Formula's weight: {metadata['formula_weight']:.2f}""")

    print(f"""------ POCKET DESCRIPTORS -----
Atom count: {pocket_descriptors['atom_count']}
Hydrophobicity: {pocket_descriptors['hydrophobic_perc']:.2f}
Number of aromatic rings: {pocket_descriptors['aromatic_count']}
Solvent accessible surface area: {pocket_descriptors['sasa']:.2f}
Gyration radius: {pocket_descriptors['gyration_radius']:.2f}""")
