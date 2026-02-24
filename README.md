# Protein Pocket Mini

Protein Pocket Mini is a small cli prototype of an application that automates the extraction and analysis of ligand binding pockets from 3D protein structures. 

## What It Does

Given one or more PDB IDs, the application executes the following pipeline:
1. **Data Retrieval:** Downloads the 3D structure (CIF format) and exact ligand metadata directly from the RCSB PDB using its GraphQL API.
2. **Pocket Identification:** Locates all ligands in the structure. It defines the "binding pocket" as all polymer atoms located within a 5.0 Å radius of the ligand.
3. **Edge Cases:**  Handle edge cases when a small molecule is unbound to the protein or when no ligands have been found.
4. **Descriptor Calculation:** Computes some geometric and chemical properties for each identified pocket:
   * **Total atom count**
   * **Hydrophobicity:** The percentage of atoms belonging to nonpolar residues.
   * **Aromaticity:** The number of atoms belonging to aromatic rings.
   * **SASA:** The total Solvent Accessible Surface Area.
   * **Gyration radius:** A measure of the pocket's spatial compactness.
5. **Reporting:** Outputs the structured analysis to the terminal, with options to export to JSON or a text file.

## Installation

Clone the repository and install the requirements.

To clone the repository run:
```bash
git clone https://github.com/melnicka/Protein-Pocket-Mini.git
```

To install the requirements run:
```bash
pip install -r requirements.txt
```

## How to Use

Run the `main.py` script from your terminal and provide one or more PDB IDs as positional arguments.

### Basic Command
```bash
python3 main.py <PDB_ID>
```

**Example:**
```bash
python3 main.py 9T1Q 9kqh
```

### Optional Arguments

You can control how the application outputs data using the following flags:

* `-t, --save-text <FILE_NAME>` : Saves the formatted terminal output to the specified text file.
* `-j, --save-json <FILE_NAME>` : Saves the complete pocket data to a JSON file.
* `-q, --quiet` : Suppresses all terminal output. Best used when saving results directly to a file.

### Usage Examples

**1. Process multiple proteins and save a copy of the output:**
Analyzes two structures, prints the results to the terminal, and appends the text to `results.txt`.
```bash
python3 main.py 9T1Q 9kqh -t results.txt
```

**2. Save machine-readable data silently:**
Analyzes a protein, saves the calculated descriptors to `data.json`, and keeps the console clean.
```bash
python3 main.py 9T1Q -qj data.json 
```
