from src.config import Config
import biotite.structure as struct
from src.data_fetching import get_protein_data 
from src.pocket import find_pockets

cfg = Config()

id="9kqh"

protein_arr, ligands = get_protein_data(id, cfg)

pocket0, pocket1 = find_pockets(protein_arr, ligands, cfg)

print(pocket0[0].shape)
print(pocket1[0].shape)
