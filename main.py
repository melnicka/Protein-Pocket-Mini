from src.config import Config
import biotite.structure as struct
from src.data_fetching import get_protein_data 
from src.pocket import find_pockets, calculate_descriptors

cfg = Config()

id="9kqh"

protein_arr, ligands = get_protein_data(id, cfg)

pockets = find_pockets(protein_arr, ligands, cfg)

print(calculate_descriptors(pockets[0]))
print('\n', calculate_descriptors(pockets[1]))

