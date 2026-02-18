from src.config import Config
import biotite.structure as struct
from src.data_fetching import get_protein_data 
cfg = Config()

id="9KQH"

at_arr, ligand_data = get_protein_data(id, cfg)

print(at_arr.shape)
print(ligand_data)
