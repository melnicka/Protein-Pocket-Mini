import os
import json
import requests
import biotite.structure.io.pdbx as pdbx
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .config import Config
    from biotite.structure import AtomArray


def get_protein_data(pdb_id: str, cfg:Config) -> tuple:
    """Acts as the main pipeline to fetch, parse, and extract protein and ligand data.

        Converts raw data into AtomArrays and dictionaries.

        Args:
            pdb_id: Unique PDB identifier.
            cfg: Config object with an attribute:
                - data_dir (str): Path to the data directory

        Returns:
            A 3-element tuple containing:
                - protein_arr (AtomArray): 3D structure of the target protein.
                - ligands (list[AtomArray]): List of extracted structures for each ligand.
                - ligand_dicts (list[dict]): List of metadata dictionaries for each ligand.
        """
    pdb_id = pdb_id.upper()
    protein_path = get_cif(pdb_id, cfg)
    cif_file = pdbx.CIFFile().read(protein_path)

    ligand_path = get_ligand_json(pdb_id, cfg)
    ligand_dicts = parse_ligand_json(ligand_path)

    protein_arr = pdbx.get_structure(cif_file, model=1)
    ligands = find_ligands(protein_arr, ligand_dicts) 

    return(protein_arr, ligands, ligand_dicts)

def get_cif(pdb_id: str, cfg: Config) -> str:
    """Downloads protein's 3D structure file in CIF format from RCSB PDB.
    Creates a unique directory for protein's data.
    If the file already exists, just returns the file path.

    Args:
        pdb_id: Unique PDB identifier.
        cfg: Config object with an attribute:
            - data_dir (str): Path to the data directory

    Returns:
        Path to the CIF file.

    Raises:
        RuntimeError: If failed to get the response from the server.
    """
    dir_path = f"{cfg.data_dir}/{pdb_id}"
    file_path = f"{dir_path}/{pdb_id}.cif"

    if os.path.exists(file_path):
        return file_path

    os.makedirs(dir_path, exist_ok=True)

    try: 
        resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif")
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to download {pdb_id}: {e}")

    with open(file_path, "wb") as f:
        f.write(resp.content)

    return file_path

def find_ligands(protein_arr: AtomArray, metadata: list[dict]) -> list[AtomArray]:
    """Finds AtomArray for each ligand present in the protein.

    Args:
        protein_arr: AtomArray of the target protein.
        metadata: List of dictionaries generated with parse_ligand_json().

    Returns:
        List of AtomArrays for each ligand.
    """
    all_ligands = []
    if len(metadata) == 0:
        return all_ligands

    for ligand in metadata:
        ligand_mask = (
            (protein_arr.chain_id == ligand['auth_asym_id']) &
            (protein_arr.res_id == int(ligand['auth_seq_id'])) &
            (protein_arr.res_name == ligand['comp_id'])
        )
        all_ligands.append(protein_arr[ligand_mask])
        
    return all_ligands

def parse_ligand_json(path: str) -> list[dict]:
    """Parses a json file generated with get_ligand_json() function.

    Args:
        path: Path to the JSON file with ligand metadata.

    Returns:
        List of dictionaries with data for each ligand. 
    """
    with open(path, 'r') as f:
        metadata = json.load(f)
        
        data = metadata['data']['entry']['nonpolymer_entities']
        ligand_list = []

        if data is None:
            return ligand_list

        for entity in data:
            comps = entity['nonpolymer_comp']['chem_comp']
            comp_dict = {
                'name': comps['name'],
                'formula': comps['formula'],
                'formula_weight': comps['formula_weight']
            }

            instances = entity["nonpolymer_entity_instances"]
            for i in range(0,len(instances)):
                ids = instances[i]["rcsb_nonpolymer_entity_instance_container_identifiers"]
                ligand_dict = {
                    "auth_asym_id": ids["auth_asym_id"],
                    "auth_seq_id": ids["auth_seq_id"],
                    "comp_id": ids["comp_id"]
                }
                ligand_dict.update(comp_dict)
                ligand_list.append(ligand_dict)

    return ligand_list

def get_ligand_json(pdb_id: str, cfg: Config) -> str:
    """Downloads ligand metadata through RCSB PDB GraphQL API and saves it in a JSON format.

    If the file already exists, just returns the file path.
    Note: This function must be ran after get_cif(). 

    Args:
        pdb_id: Unique PDB identifier.
        cfg: Config object with an attribute:
            - data_dir (str): Path to the data directory

    Returns:
        Path to the saved JSON file.

    Raises:
        RuntimeError: If failed to post the query and get the response from the server.
    """
    dir_path = f"{cfg.data_dir}/{pdb_id}"
    file_path = f"{dir_path}/{pdb_id}_ligands.json"

    if os.path.exists(file_path):
        return file_path

    query = """
{
      entry(entry_id: "%s") {
        nonpolymer_entities {
          nonpolymer_entity_instances {
            rcsb_nonpolymer_entity_instance_container_identifiers {
              auth_seq_id
              comp_id
              auth_asym_id
            }
          }
          nonpolymer_comp {
            chem_comp {
              formula_weight
              name
              formula
            }
          }
        }
      }
    }
    """ % pdb_id
    url = "https://data.rcsb.org/graphql"
    try:
        resp = requests.post(url, json={"query": query}).json()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to fetch {pdb_id} ligand metadata: {e}")

    with open(file_path, 'w') as f:
        json.dump(resp, f)

    return file_path
        
