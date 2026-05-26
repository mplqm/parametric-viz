import json
import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SVDSuperimposer import SVDSuperimposer

GEO_DIR = "geometry"
OUT_FILE = "clusters.json"
CLUSTER_THRESHOLD = 15.0  # Å — cut the dendrogram here

# ── amino-acid lookup ─────────────────────────────────────────────────────────
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H", "SEC": "U",
}

_aligner = PairwiseAligner()
_aligner.mode = "global"
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.open_gap_score = -10
_aligner.extend_gap_score = -0.5


def get_sequence(struct):
    return "".join(_AA3TO1.get(r["res_name"], "X") for r in struct["ca_trace"])


def matched_coords(struct_a, struct_b):
    """Return paired Cα arrays aligned by BLOSUM62 global sequence alignment."""
    seq_a = get_sequence(struct_a)
    seq_b = get_sequence(struct_b)
    aln = _aligner.align(seq_a, seq_b)[0]

    idx_a, idx_b = aln.indices
    mask = (idx_a != -1) & (idx_b != -1)
    ia, ib = idx_a[mask], idx_b[mask]

    ca_a = np.array([struct_a["ca_trace"][i]["coords"] for i in ia], dtype=np.float64)
    ca_b = np.array([struct_b["ca_trace"][i]["coords"] for i in ib], dtype=np.float64)
    return ca_a, ca_b, int(mask.sum())


def kabsch_fit(ca_ref, ca_mob):
    sup = SVDSuperimposer()
    sup.set(ca_ref, ca_mob)
    sup.run()
    return sup.get_rms(), sup.get_rotran()


# ── 1. Load all geometry JSON files ──────────────────────────────────────────
structures = {}
for fname in sorted(os.listdir(GEO_DIR)):
    if fname.endswith(".json"):
        with open(os.path.join(GEO_DIR, fname)) as f:
            d = json.load(f)
        structures[d["pdb_id"]] = d

pdb_ids = list(structures.keys())
n = len(pdb_ids)

# 4AKE has 2.0 Å resolution and sits outside the top-20; use first available.
ref_id = "4AKE" if "4AKE" in structures else pdb_ids[0]
print(f"Loaded {n} structures.  Reference: {ref_id}")
if ref_id != "4AKE":
    print(f"  (4AKE not in top-20 by resolution — using {ref_id} as reference)\n")

# ── 2. Align every structure to the reference (Kabsch / least-squares) ────────
print(f"Aligning all structures to {ref_id} via sequence-matched Cα…\n")
ref = structures[ref_id]
align_log = {}

for pid in pdb_ids:
    if pid == ref_id:
        align_log[pid] = {"n_aligned": ref["residue_count"], "rmsd_to_ref": 0.0}
        continue

    ca_ref, ca_mob, n_aln = matched_coords(ref, structures[pid])
    rmsd, (rot, tran) = kabsch_fit(ca_ref, ca_mob)

    all_coords = np.array(
        [r["coords"] for r in structures[pid]["ca_trace"]], dtype=np.float64
    )
    aligned = np.dot(all_coords, rot.T) + tran
    for r, c in zip(structures[pid]["ca_trace"], aligned):
        r["coords"] = [round(float(v), 4) for v in c]

    align_log[pid] = {"n_aligned": n_aln, "rmsd_to_ref": round(float(rmsd), 3)}

# ── 3. Compute 20×20 pairwise RMSD matrix ────────────────────────────────────
print("Computing pairwise RMSD matrix…")
rmsd_matrix = np.zeros((n, n))
n_aligned_matrix = np.zeros((n, n), dtype=int)

for i in range(n):
    for j in range(i + 1, n):
        ca_i, ca_j, n_aln = matched_coords(
            structures[pdb_ids[i]], structures[pdb_ids[j]]
        )
        rmsd, _ = kabsch_fit(ca_i, ca_j)
        rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd
        n_aligned_matrix[i, j] = n_aligned_matrix[j, i] = n_aln

# ── 4. Agglomerative clustering (UPGMA) ──────────────────────────────────────
condensed = squareform(rmsd_matrix)
Z = linkage(condensed, method="average")  # UPGMA
labels = fcluster(Z, t=CLUSTER_THRESHOLD, criterion="distance")

clusters: dict[int, list[str]] = {}
for pid, lbl in zip(pdb_ids, labels):
    clusters.setdefault(int(lbl), []).append(pid)

# ── 5. Save clusters.json ─────────────────────────────────────────────────────
out = {
    "reference": ref_id,
    "cluster_threshold_angstroms": CLUSTER_THRESHOLD,
    "pdb_ids": pdb_ids,
    "rmsd_matrix": [[round(v, 3) for v in row] for row in rmsd_matrix.tolist()],
    "n_clusters": len(clusters),
    "cluster_assignments": {pid: int(lbl) for pid, lbl in zip(pdb_ids, labels)},
    "clusters": {str(k): v for k, v in sorted(clusters.items())},
    "alignment_to_reference": align_log,
}
with open(OUT_FILE, "w") as f:
    json.dump(out, f, indent=2)

# ── 6. Print RMSD matrix ──────────────────────────────────────────────────────
W = 7
print(f"\nRMSD matrix (Å)  — pairwise Kabsch on BLOSUM62-matched Cα:\n")
print(f"{'':8}" + "".join(f"{pid:>{W}}" for pid in pdb_ids))
for i, pid in enumerate(pdb_ids):
    row = f"{pid:<8}" + "".join(f"{rmsd_matrix[i, j]:{W}.2f}" for j in range(n))
    print(row)

# ── Print alignment summary ───────────────────────────────────────────────────
print(f"\nAlignment to {ref_id}  (BLOSUM62 global):")
print(f"  {'ID':<8} {'Cα matched':>10}  {'RMSD to ref (Å)':>16}")
print("  " + "-" * 38)
for pid in pdb_ids:
    info = align_log[pid]
    if pid == ref_id:
        print(f"  {pid:<8} {'—':>10}  {'reference':>16}")
    else:
        print(f"  {pid:<8} {info['n_aligned']:>10}  {info['rmsd_to_ref']:>16.3f}")

# ── Print clusters ────────────────────────────────────────────────────────────
print(f"\nAgglomerative clustering (UPGMA, cut at {CLUSTER_THRESHOLD} Å):")
print(f"  → {len(clusters)} cluster(s)\n")
for lbl in sorted(clusters):
    members = clusters[lbl]
    tag = "structures" if len(members) > 1 else "structure"
    print(f"  Cluster {lbl}  ({len(members):2d} {tag}):  {', '.join(members)}")

print(f"\nSaved → {OUT_FILE}")
