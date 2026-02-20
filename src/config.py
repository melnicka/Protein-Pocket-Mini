from dataclasses import dataclass

@dataclass()
class Config:
    """Configuration settings for the protein pocket analysis pipeline.

    Attributes:
        data_dir: Path to the root directory where all unique
            protein directories will be stored. 
        radius: Search radius in Angstroms [Å] used to identify
            binding pocket atoms around a ligand's coordinates. 
    """
    data_dir: str = "data"
    radius: float = 5.0
 
