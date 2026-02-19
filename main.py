import argparse
from src.config import Config
from src.utils import display_full_info

cfg = Config()

parser = argparse.ArgumentParser()
parser.add_argument("pdb_id", help="PDB protein id.")
arg = parser.parse_args()

display_full_info(arg.pdb_id, cfg)
