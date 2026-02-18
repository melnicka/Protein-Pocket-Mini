import os
import requests
import biotite.structure as struct
import biotite.structure.io.pdbx as pdbx
from .config import Config

def get_protein_data(pdb_id, cfg:Config) -> tuple[struct.AtomArray, list[dict]]:
    path = get_cif(pdb_id, cfg)
    cif_file = pdbx.CIFFile().read(path)
    metadata = get_ligand_metadata(pdb_id)
    prot_array = pdbx.get_structure(cif_file, model=1)

    return(prot_array, metadata)

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


def get_ligand_metadata(pdb_id: str) -> list[dict]:
    query = """
{
      entry(entry_id: "%s") {
        nonpolymer_entities {
          rcsb_nonpolymer_entity_container_identifiers {
            asym_ids
          }
          nonpolymer_entity_instances {
            rcsb_nonpolymer_entity_instance_container_identifiers {
              auth_seq_id
              comp_id
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

    resp = resp['data']['entry']['nonpolymer_entities']
    ligand_list = []
    for l in resp:
        asym_ids = l["rcsb_nonpolymer_entity_container_identifiers"]["asym_ids"]
        rest = l["nonpolymer_entity_instances"]
        for i in range(0,len(asym_ids)):
            rest_containers = rest[i]["rcsb_nonpolymer_entity_instance_container_identifiers"]
            ligand_list.append({
                "asym_id": asym_ids[i],
                "auth_seq_id": rest_containers["auth_seq_id"],
                "comp_id": rest_containers["comp_id"]
            })

    return ligand_list
        
