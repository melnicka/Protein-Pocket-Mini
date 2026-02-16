from dataclasses import dataclass

@dataclass(kw_only=True)
class Config:
    data_dir: str = "data"

