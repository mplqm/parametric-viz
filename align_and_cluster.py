import json
import os
import numpy as np
from multiprocessing import Pool
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SVDSuperimposer import SVDSuperimposer

GEO_DIR           = "geometry"
FLAG_FILE         = "flagged_short_chain.json"
OUT_FILE          = "clusters.json"
CLUSTER_THRESHOLD = 15.0   # Å UPGMA cut

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H", "SEC": "U",
}

# ── Module-level globals shared with worker processes via fork ────────────────
_g_seqs   = []   # list[str]   — set before Pool creation
_g_coords = []   # list[ndarray] — set before Pool creation
_g_aln    = None
_g_sup    = None


def _worker_init():
    global _g_aln, _g_sup
    _g_aln = PairwiseAligner()
    _g_aln.mode = "global"
    _g_aln.substitution_matrix = substitution_matrices.load("BLOSUM62")
    _g_aln.open_gap_score  = -10
    _g_aln.extend_gap_score = -0.5
    _g_sup = SVDSuperimposer()


def _rmsd_worker(ij):
    i, j = ij
    aln = _g_aln.align(_g_seqs[i], _g_seqs[j])[0]
    idx_i, idx_j = aln.indices
    mask = (idx_i != -1) & (idx_j != -1)
    if mask.sum() < 3:
        return 99.9
    ci = _g_coords[i][idx_i[mask]]
    cj = _g_coords[j][idx_j[mask]]
    _g_sup.set(ci, cj)
    _g_sup.run()
    return float(_g_sup.get_rms())


def get_sequence(struct):
    return "".join(_AA3TO1.get(r["res_name"], "X") for r in struct["ca_trace"])


def kabsch_align(struct, ref):
    seq_r = get_sequence(ref)
    seq_m = get_sequence(struct)
    aln   = _aligner_main.align(seq_r, seq_m)[0]
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

    return int(mask.sum()), float(sup.get_rms())


_aligner_main = PairwiseAligner()
_aligner_main.mode = "global"
_aligner_main.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner_main.open_gap_score  = -10
_aligner_main.extend_gap_score = -0.5


if __name__ == "__main__":
    # ── 1. Load structures, skip short-chain flagged ──────────────────────────
    short_chain_ids = set()
    if os.path.exists(FLAG_FILE):
        with open(FLAG_FILE) as f:
            flag_data = json.load(f)
        short_chain_ids = set(flag_data.get("flagged", {}).keys())
        print(f"Short-chain flags loaded: {len(short_chain_ids)} structures will be skipped")

    structures = {}
    for fname in sorted(os.listdir(GEO_DIR)):
        if not fname.endswith(".json"):
            continue
        pid = fname[:-5]
        if pid in short_chain_ids:
            continue
        with open(os.path.join(GEO_DIR, fname)) as f:
            structures[pid] = json.load(f)

    pdb_ids = list(structures.keys())
    N       = len(pdb_ids)
    ref_id  = "4AKE" if "4AKE" in structures else pdb_ids[0]
    print(f"Loaded {N} structures (after short-chain exclusion).  Reference: {ref_id}\n")

    # ── 2. Align all structures to reference (sequential) ────────────────────
    print(f"Aligning {N} structures to {ref_id}…")
    ref = structures[ref_id]
    align_log = {ref_id: {"n_aligned": ref["residue_count"], "rmsd_to_ref": 0.0}}

    for i, pid in enumerate(pdb_ids):
        if pid == ref_id:
            continue
        n_aln, rmsd = kabsch_align(structures[pid], ref)
        align_log[pid] = {"n_aligned": n_aln, "rmsd_to_ref": round(rmsd, 3)}
        if (i + 1) % 50 == 0:
            print(f"  aligned {i + 1}/{N}")

    print(f"  aligned {N}/{N}  ✓")

    # ── 3. Flag RMSD outliers (> 3× median rmsd_to_ref) ──────────────────────
    rmsd_to_ref = np.array(
        [align_log[p]["rmsd_to_ref"] for p in pdb_ids if p != ref_id]
    )
    median_rmsd     = float(np.median(rmsd_to_ref))
    outlier_cutoff  = 3.0 * median_rmsd
    rmsd_outliers   = {
        pid for pid in pdb_ids
        if pid != ref_id and align_log[pid]["rmsd_to_ref"] > outlier_cutoff
    }
    print(f"\nRMSD-to-reference: median={median_rmsd:.2f} Å  outlier cutoff (3×)={outlier_cutoff:.2f} Å")
    print(f"  RMSD outliers: {len(rmsd_outliers)}")
    if rmsd_outliers:
        for pid in sorted(rmsd_outliers):
            print(f"    {pid:<8} rmsd={align_log[pid]['rmsd_to_ref']:.2f} Å")

    active_ids = [p for p in pdb_ids if p not in rmsd_outliers]
    A          = len(active_ids)
    n_pairs    = A * (A - 1) // 2
    n_cpu      = min(8, os.cpu_count() or 1)
    print(f"\nActive structures for clustering: {A}  ({n_pairs} pairs)  using {n_cpu} cores")

    # ── 4. Parallel pairwise RMSD ─────────────────────────────────────────────
    print("Computing pairwise RMSD matrix…")
    _g_seqs[:]   = [get_sequence(structures[p]) for p in active_ids]
    _g_coords[:] = [
        np.array([r["coords"] for r in structures[p]["ca_trace"]], dtype=np.float64)
        for p in active_ids
    ]

    pairs      = [(i, j) for i in range(A) for j in range(i + 1, A)]
    chunksize  = max(1, len(pairs) // (n_cpu * 8))

    with Pool(processes=n_cpu, initializer=_worker_init) as pool:
        flat = pool.map(_rmsd_worker, pairs, chunksize=chunksize)

    rmsd_matrix = np.zeros((A, A))
    for (i, j), v in zip(pairs, flat):
        rmsd_matrix[i, j] = rmsd_matrix[j, i] = v

    print(f"  done  — min={rmsd_matrix[rmsd_matrix > 0].min():.2f}  "
          f"max={rmsd_matrix.max():.2f}  mean={rmsd_matrix[rmsd_matrix > 0].mean():.2f} Å")

    # ── 5. UPGMA clustering ───────────────────────────────────────────────────
    condensed = squareform(rmsd_matrix)
    Z         = linkage(condensed, method="average")
    labels    = fcluster(Z, t=CLUSTER_THRESHOLD, criterion="distance")

    clusters: dict[int, list[str]] = {}
    for pid, lbl in zip(active_ids, labels):
        clusters.setdefault(int(lbl), []).append(pid)

    print(f"\nClustering (UPGMA, {CLUSTER_THRESHOLD} Å):  {len(clusters)} clusters")
    for lbl in sorted(clusters):
        print(f"  C{lbl:<3}  {len(clusters[lbl]):>3} members:  {', '.join(clusters[lbl][:8])}"
              + ("…" if len(clusters[lbl]) > 8 else ""))

    # ── 6. Save clusters.json ─────────────────────────────────────────────────
    condensed_list = [[round(v, 3) for v in row] for row in rmsd_matrix.tolist()]
    out = {
        "reference":                  ref_id,
        "cluster_threshold_angstroms": CLUSTER_THRESHOLD,
        "n_structures":               A,
        "n_clusters":                 len(clusters),
        "pdb_ids":                    active_ids,
        "rmsd_matrix":                condensed_list,
        "cluster_assignments":        {pid: int(lbl) for pid, lbl in zip(active_ids, labels)},
        "clusters":                   {str(k): v for k, v in sorted(clusters.items())},
        "alignment_to_reference":     align_log,
        "excluded": {
            "short_chain": sorted(short_chain_ids),
            "rmsd_outliers": sorted(rmsd_outliers),
        },
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → {OUT_FILE}")
