from src.config import Config
from src.utils import get_cif

cfg = Config()

pdb_path = get_cif("2PGH", cfg)
print(pdb_path)
