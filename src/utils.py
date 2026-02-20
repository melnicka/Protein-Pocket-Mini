from .data_fetching import get_protein_data
from .pocket import find_pockets, calculate_descriptors
from .config import Config

def display_full_info(pdb_id: str, cfg: Config):
    pdb_id = pdb_id.upper()
    protein_arr, ligands, ligand_dict = get_protein_data(pdb_id, cfg)
    num_pockets = len(ligands)
    pockets = find_pockets(protein_arr, ligands, cfg)

    print(f"\n{'=' * 60}")
    print(f" PROTEIN: {pdb_id} | TOTAL POCKETS FOUND: {num_pockets}")
    print(f"{'=' * 60}") 

    for i in range(0, num_pockets):
        print(f"\n----------------[ POCKET {i+1} ]----------------")
        pocket = pockets[i]
        display_ligand_info(ligand_dict[i])

        if pocket is None:
            print("[ Ligand not bount to the protein ]")
        else:
            descriptors = calculate_descriptors(pocket)
            display_pocket_info(descriptors)
    print("")

def display_ligand_info(ligand_dict: dict):
    print(f"""[ Ligand Info ]
    • Name: {ligand_dict['name']}
    • Formula: {ligand_dict['formula']}
    • Formula's weight: {ligand_dict['formula_weight']:.2f}
    """)

def display_pocket_info(pocket_descriptors: dict):
    print(f"""[ Pocket Descriptors ]
    • Atom count: {pocket_descriptors['atom_count']}
    • Hydrophobicity: {pocket_descriptors['hydrophobic_perc']:.2f}
    • Number of aromatic rings: {pocket_descriptors['aromatic_count']}
    • Solvent accessible surface area: {pocket_descriptors['sasa']:.2f}
    • Gyration radius: {pocket_descriptors['gyration_radius']:.2f}""")
