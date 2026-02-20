import argparse
from src.config import Config
from src.utils import display_full_info

if __name__ == "__main__":

    cfg = Config()

    parser = argparse.ArgumentParser(
        description= "Find all ligand binding pockets and calculate its descriptors."
    )
    parser.add_argument(
        "pdb_ids",
        nargs="+",
        type=str,
        help="List of PDB IDs to process (e.g., 9T1Q 9kqh 9ssd)"
    )
    args = parser.parse_args()

    for pdb_id in args.pdb_ids:
        display_full_info(pdb_id, cfg)
