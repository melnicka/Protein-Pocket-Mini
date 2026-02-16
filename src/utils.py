import os
import requests
from .config import Config

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


