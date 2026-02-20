from .data_fetching import get_protein_data
from .pocket import find_pockets, calculate_descriptors
import argparse
from .config import Config

def get_full_info(pdb_id: str, cfg: Config) -> dict:
    pdb_id = pdb_id.upper()
    protein_arr, ligands, ligand_dicts = get_protein_data(pdb_id, cfg)
    num_pockets = len(ligands)
    pockets = find_pockets(protein_arr, ligands, cfg)
    full_data = {'pdb_id': pdb_id, 'num_pockets': num_pockets}

    for i in range(0, num_pockets):
        full_data[f'pocket{i+1}'] = {}
        full_data[f'pocket{i+1}']['ligand_info'] = ligand_dicts[i]
        if pockets[i] is None:
            full_data[f'pocket{i+1}']['descriptors'] = None
        else:
            descriptors = calculate_descriptors(pockets[i])
            full_data[f'pocket{i+1}']['descriptors'] = descriptors

    return full_data

def display_full_info(full_data: dict):

    print(f"\n{'=' * 60}")
    print(f" PROTEIN: {full_data['pdb_id']} | TOTAL POCKETS FOUND: {full_data['num_pockets']}")
    print(f"{'=' * 60}") 

    for i in range(0, full_data['num_pockets']):
        print(f"\n----------------[ POCKET {i+1} ]----------------")
        display_ligand_info(full_data[f'pocket{i+1}']['ligand_info'])

        if full_data[f'pocket{i+1}']['descriptors'] is None:
            print("[ Ligand not bount to the protein ]")
        else:
            display_pocket_info(full_data[f'pocket{i+1}']['descriptors'] )
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

def parse_args():
    parser = argparse.ArgumentParser(
        description="""Find all ligand binding pockets - get information about
ligands and calculate descriptors for each pocket. If the ligand is not bound to the
protein, i.e., surrounded by solvent alone, only the information about the ligand
will be displayed.""")

    parser.add_argument(
        "pdb_ids",
        nargs="+",
        type=str,
        help="List of PDB IDs to process (e.g., 9T1Q 9kqh 9ssd)"
    )

    parser.add_argument(
        "--save-text", "-t",
        metavar='FILE_NAME',
        type=str,
        help="Optional: Saves terminal output to a specified text file."
    )

    parser.add_argument(
        "--save-json", "-j",
        type=str,
        metavar='FILE_NAME',
        help="Optional: Saves full protein data to a specified json file."
    )

    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Optional: Disable terminal output."
    )

    args = parser.parse_args()
    
    return args


