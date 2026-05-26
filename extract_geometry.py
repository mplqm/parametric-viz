import json
import os
import numpy as np
from Bio import PDB

PDB_DIR = "pdb_files"
OUT_DIR = "geometry"
os.makedirs(OUT_DIR, exist_ok=True)

parser = PDB.PDBParser(QUIET=True)

print(f"{'PDB ID':<8} {'Chains':<8} {'Residues':>8}")
print("-" * 26)

for fname in sorted(os.listdir(PDB_DIR)):
    if not fname.endswith(".pdb"):
        continue

    pdb_id = fname[:-4]
    structure = parser.get_structure(pdb_id, os.path.join(PDB_DIR, fname))

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  # skip HETATM and waters
                    continue
                if "CA" not in residue:
                    continue
                ca = residue["CA"]
                residues.append({
                    "chain": chain.id,
                    "res_seq": residue.id[1],
                    "res_name": residue.resname,
                    "coords": ca.get_vector().get_array().tolist(),
                })
        break  # first model only

    if not residues:
        print(f"{pdb_id:<8} {'—':<8} {'0':>8}  (no CA atoms found)")
        continue

    coords = np.array([r["coords"] for r in residues])
    centroid = coords.mean(axis=0)
    centered = (coords - centroid).tolist()

    for r, c in zip(residues, centered):
        r["coords"] = [round(v, 4) for v in c]

    chains = sorted({r["chain"] for r in residues})
    out = {
        "pdb_id": pdb_id,
        "chains": chains,
        "centroid_original": [round(v, 4) for v in centroid.tolist()],
        "residue_count": len(residues),
        "ca_trace": residues,
    }

    with open(os.path.join(OUT_DIR, f"{pdb_id}.json"), "w") as f:
        json.dump(out, f, indent=2)

    print(f"{pdb_id:<8} {','.join(chains):<8} {len(residues):>8}")
