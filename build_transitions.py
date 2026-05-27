import json
import os
import shutil
import time
import numpy as np
from itertools import combinations
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SVDSuperimposer import SVDSuperimposer

GEO_DIR          = "geometry"
CLUSTERS_FILE    = "clusters.json"
OUT_DIR          = "transitions"
TMP_DIR          = "transitions_tmp"
N_FRAMES         = 60
MIN_RESIDUES     = 20
MIN_CLUSTER_SIZE = 2   # skip singleton↔singleton transitions

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
_aligner.open_gap_score   = -10
_aligner.extend_gap_score = -0.5


def get_sequence(struct):
    return "".join(_AA3TO1.get(r["res_name"], "X") for r in struct["ca_trace"])


def kabsch_align(struct, ref):
    seq_r = get_sequence(ref)
    seq_m = get_sequence(struct)
    aln   = _aligner.align(seq_r, seq_m)[0]
    ir, im = aln.indices
    mask   = (ir != -1) & (im != -1)

    ca_r = np.array([ref["ca_trace"][i]["coords"]    for i in ir[mask]], dtype=np.float64)
    ca_m = np.array([struct["ca_trace"][i]["coords"] for i in im[mask]], dtype=np.float64)

    sup = SVDSuperimposer()
    sup.set(ca_r, ca_m)
    sup.run()
    rot, tran = sup.get_rotran()

    all_c   = np.array([r["coords"] for r in struct["ca_trace"]], dtype=np.float64)
    aligned = np.dot(all_c, rot.T) + tran
    for r, c in zip(struct["ca_trace"], aligned):
        r["coords"] = [round(float(v), 4) for v in c]

    return float(sup.get_rms())


def ref_mapping(struct, ref):
    seq_r = get_sequence(ref)
    seq_m = get_sequence(struct)
    aln   = _aligner.align(seq_r, seq_m)[0]
    ir, im = aln.indices
    mask   = (ir != -1) & (im != -1)
    return {
        int(r): np.array(struct["ca_trace"][m]["coords"], dtype=np.float64)
        for r, m in zip(ir[mask], im[mask])
    }


def smoothstep(t):
    return t * t * (3.0 - 2.0 * t)


t_start = time.time()

# ── 1. Load ───────────────────────────────────────────────────────────────────
structures = {}
for fname in sorted(os.listdir(GEO_DIR)):
    if fname.endswith(".json"):
        with open(os.path.join(GEO_DIR, fname)) as f:
            d = json.load(f)
        structures[d["pdb_id"]] = d

with open(CLUSTERS_FILE) as f:
    cluster_data = json.load(f)

ref_id   = cluster_data["reference"]
ref      = structures[ref_id]
clusters = cluster_data["clusters"]
print(f"Reference: {ref_id}   Clusters: {len(clusters)}   Frames/transition: {N_FRAMES}\n")

# ── 2. Align all to reference ─────────────────────────────────────────────────
active_pids = {pid for members in clusters.values() for pid in members}
print(f"Aligning {len(active_pids)} structures to {ref_id}…")
for i, pid in enumerate(sorted(active_pids)):
    if pid != ref_id:
        kabsch_align(structures[pid], ref)
    if (i + 1) % 100 == 0:
        print(f"  {i + 1}/{len(active_pids)}")
print(f"  {len(active_pids)}/{len(active_pids)}  ✓")

t_align = time.time()

ref_self_map = {i: np.array(r["coords"], dtype=np.float64)
                for i, r in enumerate(ref["ca_trace"])}

# ── 3. Cluster centroids ──────────────────────────────────────────────────────
print("\nComputing cluster centroids…")
centroids = {}
for label, members in clusters.items():
    maps = []
    for pid in members:
        if pid == ref_id:
            maps.append(ref_self_map)
        else:
            maps.append(ref_mapping(structures[pid], ref))

    common   = sorted(set.intersection(*[set(m.keys()) for m in maps]))
    centroid = {pos: np.mean([m[pos] for m in maps], axis=0) for pos in common}
    centroids[label] = centroid

    res_ids = [ref["ca_trace"][p]["res_seq"] for p in common]
    span    = f"{res_ids[0]}–{res_ids[-1]}" if res_ids else "—"
    print(f"  C{label:<3} {len(members):>3} members  {len(common):>4} common Cα  ref {span}")

t_centroids = time.time()

# ── 4. Generate transitions → write to TMP_DIR first (atomic swap) ───────────
n_pairs = len(clusters) * (len(clusters) - 1) // 2
print(f"\nGenerating up to {n_pairs:,} pairwise transitions…")

if os.path.exists(TMP_DIR):
    shutil.rmtree(TMP_DIR)
os.makedirs(TMP_DIR)

cluster_labels = sorted(clusters.keys(), key=int)
index_entries  = []
skipped        = 0
written        = 0

for lbl_a, lbl_b in combinations(cluster_labels, 2):
    if len(clusters[lbl_a]) < MIN_CLUSTER_SIZE and len(clusters[lbl_b]) < MIN_CLUSTER_SIZE:
        skipped += 1
        continue

    cent_a = centroids[lbl_a]
    cent_b = centroids[lbl_b]
    common = sorted(set(cent_a) & set(cent_b))
    if len(common) < MIN_RESIDUES:
        skipped += 1
        continue

    res_seq_list = [ref["ca_trace"][p]["res_seq"] for p in common]
    coords_a     = np.array([cent_a[p] for p in common])
    coords_b     = np.array([cent_b[p] for p in common])

    frames = []
    for k in range(N_FRAMES):
        t      = smoothstep(k / (N_FRAMES - 1))
        interp = (1.0 - t) * coords_a + t * coords_b
        frames.append([[round(float(v), 3) for v in xyz] for xyz in interp])

    tid = f"C{lbl_a}-C{lbl_b}"
    transition = {
        "from_cluster": lbl_a,
        "to_cluster":   lbl_b,
        "from_members": clusters[lbl_a],
        "to_members":   clusters[lbl_b],
        "n_residues":   len(common),
        "res_seq":      res_seq_list,
        "frames":       frames,
    }

    with open(os.path.join(TMP_DIR, f"{tid}.json"), "w") as f:
        json.dump(transition, f, separators=(",", ":"))

    index_entries.append({
        "id":           tid,
        "from_cluster": lbl_a,
        "to_cluster":   lbl_b,
        "from_members": clusters[lbl_a],
        "to_members":   clusters[lbl_b],
        "n_residues":   len(common),
        "res_seq":      res_seq_list,
    })
    written += 1
    if written % 1000 == 0:
        print(f"  {written:,} transitions written…")

# ── 5. Write index.json into TMP_DIR ─────────────────────────────────────────
centroid_summary = {}
for label, cent in centroids.items():
    pos_sorted = sorted(cent.keys())
    centroid_summary[label] = {
        "members":    clusters[label],
        "n_residues": len(pos_sorted),
        "res_seq":    [ref["ca_trace"][p]["res_seq"] for p in pos_sorted],
        "coords":     [[round(float(v), 3) for v in cent[p]] for p in pos_sorted],
    }

index = {
    "reference":     ref_id,
    "n_clusters":    len(clusters),
    "n_frames":      N_FRAMES,
    "interpolation": "smoothstep",
    "centroids":     centroid_summary,
    "transitions":   index_entries,
}
with open(os.path.join(TMP_DIR, "index.json"), "w") as f:
    json.dump(index, f, separators=(",", ":"))

t_transitions = time.time()

# ── 6. Atomic swap: replace transitions/ with transitions_tmp/ ────────────────
if os.path.exists(OUT_DIR):
    shutil.rmtree(OUT_DIR)
os.rename(TMP_DIR, OUT_DIR)
print(f"\nAtomically replaced {OUT_DIR}/")

# ── 7. Summary ────────────────────────────────────────────────────────────────
idx_kb   = os.path.getsize(os.path.join(OUT_DIR, "index.json")) / 1024
total_mb = sum(
    os.path.getsize(os.path.join(OUT_DIR, f))
    for f in os.listdir(OUT_DIR)
) / 1024 / 1024

t_end = time.time()

print(f"\n{'═'*54}")
print(f"  Transitions written:   {written:>6,}")
print(f"  Skipped:               {skipped:>6,}")
print(f"  index.json:            {idx_kb:>6.0f} KB")
print(f"  Total transitions/:    {total_mb:>6.0f} MB  ({len(os.listdir(OUT_DIR))} files)")
print(f"\n  Align:                 {t_align      - t_start:.1f}s")
print(f"  Centroids:             {t_centroids  - t_align:.1f}s")
print(f"  Transitions:           {t_transitions- t_centroids:.1f}s")
print(f"  Total:                 {t_end        - t_start:.1f}s")
print(f"{'═'*54}")

# ── 8. Summary table (first 25) ───────────────────────────────────────────────
print(f"\n{'Transition':<14} {'Residues':>9}  {'Centroid RMSD (Å)':>18}")
print("─" * 46)
for entry in index_entries[:25]:
    la, lb = entry["from_cluster"], entry["to_cluster"]
    common = sorted(set(centroids[la]) & set(centroids[lb]))
    ca_arr = np.array([centroids[la][p] for p in common])
    cb_arr = np.array([centroids[lb][p] for p in common])
    rmsd   = float(np.sqrt(np.mean(np.sum((ca_arr - cb_arr) ** 2, axis=1))))
    print(f"  {entry['id']:<12} {entry['n_residues']:>9}  {rmsd:>18.2f}")
if len(index_entries) > 25:
    print(f"  … {len(index_entries) - 25:,} more transitions")
