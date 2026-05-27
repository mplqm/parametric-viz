import json
import os
import numpy as np
from Bio import PDB

PDB_DIR  = "pdb_files"
OUT_DIR  = "geometry"
FLAG_FILE = "flagged_short_chain.json"
os.makedirs(OUT_DIR, exist_ok=True)

existing_geo = {f[:-5] for f in os.listdir(OUT_DIR) if f.endswith(".json")}
parser = PDB.PDBParser(QUIET=True)

new_count = skipped = error_count = 0

for fname in sorted(os.listdir(PDB_DIR)):
    if not fname.endswith(".pdb"):
        continue
    pdb_id = fname[:-4]

    if pdb_id in existing_geo:
        skipped += 1
        continue

    structure = parser.get_structure(pdb_id, os.path.join(PDB_DIR, fname))
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CA" not in residue:
                    continue
                residues.append({
                    "chain":    chain.id,
                    "res_seq":  residue.id[1],
                    "res_name": residue.resname,
                    "coords":   residue["CA"].get_vector().get_array().tolist(),
                })
        break  # first model only

    if not residues:
        print(f"  {pdb_id:<8} SKIP — no CA atoms")
        error_count += 1
        continue

    coords   = np.array([r["coords"] for r in residues])
    centroid = coords.mean(axis=0)
    centered = (coords - centroid).tolist()
    for r, c in zip(residues, centered):
        r["coords"] = [round(v, 4) for v in c]

    chains = sorted({r["chain"] for r in residues})
    out = {
        "pdb_id":            pdb_id,
        "chains":            chains,
        "centroid_original": [round(v, 4) for v in centroid.tolist()],
        "residue_count":     len(residues),
        "ca_trace":          residues,
    }
    with open(os.path.join(OUT_DIR, f"{pdb_id}.json"), "w") as f:
        json.dump(out, f, indent=2)

    new_count += 1
    print(f"  {pdb_id:<8} chains={','.join(chains):<6} residues={len(residues):>4}")

print(f"\n  New: {new_count}   Skipped (existing): {skipped}   Errors: {error_count}")

# ── Chain A length flagging ───────────────────────────────────────────────────
print("\nComputing chain A length statistics…")

chain_a_lengths = {}
for fname in sorted(os.listdir(OUT_DIR)):
    if not fname.endswith(".json"):
        continue
    with open(os.path.join(OUT_DIR, fname)) as f:
        d = json.load(f)
    n_chain_a = sum(1 for r in d["ca_trace"] if r["chain"] == "A")
    if n_chain_a > 0:
        chain_a_lengths[d["pdb_id"]] = n_chain_a

if not chain_a_lengths:
    print("  No structures have chain A — skipping chain-A flagging.")
    flagged = {}
else:
    lengths   = list(chain_a_lengths.values())
    median_ca = float(np.median(lengths))
    threshold = median_ca * 0.5

    flagged = {
        pid: n
        for pid, n in chain_a_lengths.items()
        if n < threshold
    }

    print(f"  Structures with chain A: {len(chain_a_lengths)}")
    print(f"  Median chain A length:   {median_ca:.0f} residues")
    print(f"  Threshold (50%%):         {threshold:.0f} residues")
    print(f"  Flagged as short:        {len(flagged)}")
    if flagged:
        for pid, n in sorted(flagged.items(), key=lambda x: x[1]):
            print(f"    {pid:<8} chain A = {n} residues")

    with open(FLAG_FILE, "w") as f:
        json.dump({
            "median_chain_a_length": round(median_ca, 1),
            "threshold_length":      round(threshold, 1),
            "flagged":               {pid: n for pid, n in sorted(flagged.items())},
        }, f, indent=2)
    print(f"\nSaved → {FLAG_FILE}")
