import biotite.structure as struct
from biotite.structure import AtomArray, CellList
from .config import Config



def find_pockets(
    protein_arr: AtomArray,
    ligand_list: list[AtomArray],
    cfg: Config
) -> list[AtomArray]:
    cell_list = CellList(protein_arr, cell_size=cfg.radius)
    pockets = []
    for ligand_arr in ligand_list:
        coords = ligand_arr.coord
        raw_indicies = cell_list.get_atoms(coords, radius=cfg.radius)
        indicies = raw_indicies[raw_indicies != -1]

        raw_pocket = protein_arr[indicies]
        pocket_mask = (
            (struct.filter_polymer(raw_pocket)) |
            (struct.filter_polymer(raw_pocket, pol_type='nucleotide')) |
            (struct.filter_polymer(raw_pocket, pol_type='carbohydrate'))
        )
        pocket = raw_pocket[pocket_mask]

        pockets.append((pocket, ligand_arr))

    return pockets


