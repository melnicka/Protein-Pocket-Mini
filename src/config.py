from dataclasses import dataclass

@dataclass()
class Config:
    data_dir: str = "data"
    radius: float = 5.0
 
