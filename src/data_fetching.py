import os
import requests
from biotite.structure import AtomArray
import biotite.structure.io.pdbx as pdbx
from .config import Config

def get_protein_data(pdb_id, cfg:Config) -> tuple[AtomArray, list[AtomArray]]:
    path = get_cif(pdb_id, cfg)
    cif_file = pdbx.CIFFile().read(path)
    metadata = get_ligand_metadata(pdb_id)
    protein_arr = pdbx.get_structure(cif_file, model=1)
    ligands = find_ligands(protein_arr, metadata)

    return(protein_arr, ligands)

def get_cif(pdb_id: str, cfg: Config) -> str:
    pdb_id = pdb_id.upper()
    file_path = f"{cfg.data_dir}/{pdb_id}.cif"
    if os.path.exists(file_path):
        return file_path

    if not os.path.exists(cfg.data_dir):
        os.makedirs(cfg.data_dir)

    try: 
        resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif")
        resp.raise_for_status()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to download {pdb_id}: {e}")

    with open(file_path, "wb") as f:
        f.write(resp.content)

    return file_path

def find_ligands(protein_arr: AtomArray, metadata: list[dict]):
    all_ligands = []
    for ligand in metadata:
        ligand_mask = (
            (protein_arr.chain_id == ligand['auth_asym_id']) &
            (protein_arr.res_id == int(ligand['auth_seq_id'])) &
            (protein_arr.res_name == ligand['comp_id'])
        )
        all_ligands.append(protein_arr[ligand_mask])
        
    return all_ligands

def get_ligand_metadata(pdb_id: str) -> list[dict]:
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
        }
      }
    }
    """ % pdb_id
    url = "https://data.rcsb.org/graphql"
    try:
        resp = requests.post(url, json={"query": query}).json()
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to fetch {pdb_id} ligand metadata: {e}")

    data = resp['data']['entry']['nonpolymer_entities']
    ligand_list = []
    for entity in data:
        instances = entity["nonpolymer_entity_instances"]
        for i in range(0,len(instances)):
            ids = instances[i]["rcsb_nonpolymer_entity_instance_container_identifiers"]
            ligand_list.append({
                "auth_asym_id": ids["auth_asym_id"],
                "auth_seq_id": ids["auth_seq_id"],
                "comp_id": ids["comp_id"]
            })

    return ligand_list
        
