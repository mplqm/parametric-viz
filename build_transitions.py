import json
import os
import numpy as np
from itertools import combinations
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SVDSuperimposer import SVDSuperimposer

GEO_DIR = "geometry"
CLUSTERS_FILE = "clusters.json"
OUT_FILE = "transitions.json"
N_FRAMES = 60

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


def kabsch_align(struct, ref):
    """Align struct to ref in-place using BLOSUM62-matched Cα. Returns RMSD."""
    seq_ref = get_sequence(ref)
    seq_mob = get_sequence(struct)
    aln = _aligner.align(seq_ref, seq_mob)[0]
    idx_ref, idx_mob = aln.indices
    mask = (idx_ref != -1) & (idx_mob != -1)

    ca_ref = np.array([ref["ca_trace"][i]["coords"] for i in idx_ref[mask]], dtype=np.float64)
    ca_mob = np.array([struct["ca_trace"][i]["coords"] for i in idx_mob[mask]], dtype=np.float64)

    sup = SVDSuperimposer()
    sup.set(ca_ref, ca_mob)
    sup.run()
    rot, tran = sup.get_rotran()

    all_coords = np.array([r["coords"] for r in struct["ca_trace"]], dtype=np.float64)
    aligned = np.dot(all_coords, rot.T) + tran
    for r, c in zip(struct["ca_trace"], aligned):
        r["coords"] = [round(float(v), 4) for v in c]

    return sup.get_rms()


def ref_mapping(struct, ref):
    """Return dict: ref_array_index → aligned Cα coords in ref frame."""
    seq_ref = get_sequence(ref)
    seq_mob = get_sequence(struct)
    aln = _aligner.align(seq_ref, seq_mob)[0]
    idx_ref, idx_mob = aln.indices
    mask = (idx_ref != -1) & (idx_mob != -1)

    return {
        int(ir): np.array(struct["ca_trace"][im]["coords"], dtype=np.float64)
        for ir, im in zip(idx_ref[mask], idx_mob[mask])
    }


def smoothstep(t: float) -> float:
    """Cubic ease-in/ease-out: 0→0, 1→1, zero derivative at endpoints."""
    return t * t * (3.0 - 2.0 * t)


# ── 1. Load structures and clusters ──────────────────────────────────────────
structures = {}
for fname in sorted(os.listdir(GEO_DIR)):
    if fname.endswith(".json"):
        with open(os.path.join(GEO_DIR, fname)) as f:
            d = json.load(f)
        structures[d["pdb_id"]] = d

with open(CLUSTERS_FILE) as f:
    cluster_data = json.load(f)

ref_id = cluster_data["reference"]
ref = structures[ref_id]
clusters = cluster_data["clusters"]  # {"1": ["3QMA", "8ESU", ...], ...}
print(f"Reference: {ref_id}   Clusters: {len(clusters)}   Frames/transition: {N_FRAMES}\n")

# ── 2. Align all structures to reference ──────────────────────────────────────
print("Aligning all structures to reference…")
for pid, struct in structures.items():
    if pid != ref_id:
        kabsch_align(struct, ref)

# Reference maps to itself (identity)
ref_self_map = {i: np.array(r["coords"], dtype=np.float64)
                for i, r in enumerate(ref["ca_trace"])}

# ── 3. Compute cluster centroids ──────────────────────────────────────────────
print("Computing cluster centroids…\n")

centroids = {}  # cluster_label -> {ref_idx: avg_coords}

for label, members in clusters.items():
    maps = []
    for pid in members:
        if pid == ref_id:
            maps.append(ref_self_map)
        else:
            maps.append(ref_mapping(structures[pid], ref))

    # Positions present in every member of this cluster
    common_pos = sorted(set.intersection(*[set(m.keys()) for m in maps]))

    centroid = {}
    for pos in common_pos:
        centroid[pos] = np.mean([m[pos] for m in maps], axis=0)

    centroids[label] = centroid

    res_ids = [ref["ca_trace"][p]["res_seq"] for p in common_pos]
    span = f"{res_ids[0]}–{res_ids[-1]}" if res_ids else "—"
    print(f"  Cluster {label:>2}  members: {', '.join(members):<35}  "
          f"centroid residues: {len(common_pos):>3}  ref span: {span}")

# ── 4. Generate pairwise transitions ─────────────────────────────────────────
print(f"\nGenerating {len(clusters)*(len(clusters)-1)//2} pairwise transitions "
      f"({N_FRAMES} frames each)…")

cluster_labels = sorted(clusters.keys(), key=int)
transitions = []

for lbl_a, lbl_b in combinations(cluster_labels, 2):
    cent_a = centroids[lbl_a]
    cent_b = centroids[lbl_b]

    common_pos = sorted(set(cent_a) & set(cent_b))
    if not common_pos:
        print(f"  WARNING: clusters {lbl_a}↔{lbl_b} share no common residues — skipped")
        continue

    res_seq_list = [ref["ca_trace"][p]["res_seq"] for p in common_pos]
    coords_a = np.array([cent_a[p] for p in common_pos])  # (N, 3)
    coords_b = np.array([cent_b[p] for p in common_pos])  # (N, 3)

    frames = []
    for k in range(N_FRAMES):
        t = smoothstep(k / (N_FRAMES - 1))
        interp = (1.0 - t) * coords_a + t * coords_b
        frames.append([[round(float(v), 3) for v in xyz] for xyz in interp])

    transitions.append({
        "from_cluster": lbl_a,
        "to_cluster": lbl_b,
        "from_members": clusters[lbl_a],
        "to_members": clusters[lbl_b],
        "n_residues": len(common_pos),
        "res_seq": res_seq_list,
        "frames": frames,
    })

# ── 5. Save transitions.json ──────────────────────────────────────────────────
centroid_summary = {}
for label, cent in centroids.items():
    common_pos = sorted(cent.keys())
    centroid_summary[label] = {
        "members": clusters[label],
        "n_residues": len(common_pos),
        "res_seq": [ref["ca_trace"][p]["res_seq"] for p in common_pos],
        "coords": [[round(float(v), 3) for v in cent[p]] for p in common_pos],
    }

out = {
    "reference": ref_id,
    "n_clusters": len(clusters),
    "n_frames": N_FRAMES,
    "interpolation": "smoothstep",
    "centroids": centroid_summary,
    "transitions": transitions,
}

with open(OUT_FILE, "w") as f:
    json.dump(out, f, separators=(",", ":"))  # compact — no whitespace

size_kb = os.path.getsize(OUT_FILE) / 1024
print(f"\nSaved → {OUT_FILE}  ({size_kb:.0f} KB,  {len(transitions)} transitions)")

# ── 6. Summary table ──────────────────────────────────────────────────────────
print(f"\n{'Transition':<14} {'Residues':>9}  {'RMSD A↔B (Å)':>14}")
print("-" * 42)
for t in transitions:
    la, lb = t["from_cluster"], t["to_cluster"]
    common = sorted(set(centroids[la]) & set(centroids[lb]))
    ca = np.array([centroids[la][p] for p in common])
    cb = np.array([centroids[lb][p] for p in common])
    # centroid-to-centroid RMSD (no fitting — both in the reference frame)
    rmsd = float(np.sqrt(np.mean(np.sum((ca - cb) ** 2, axis=1))))
    tag = f"C{la}→C{lb}"
    print(f"  {tag:<12} {t['n_residues']:>9}  {rmsd:>14.2f}")
