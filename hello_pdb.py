import urllib.request
from Bio import PDB

PDB_ID = "4AKE"
pdb_file = f"{PDB_ID}.pdb"

url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
urllib.request.urlretrieve(url, pdb_file)
print(f"Downloaded {pdb_file}")

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure(PDB_ID, pdb_file)

atom_count = sum(1 for _ in structure.get_atoms())
print(f"Total atoms in {PDB_ID}: {atom_count}")
