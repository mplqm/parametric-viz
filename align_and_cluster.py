import json
import os
import time
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
N_LANDMARKS       = 20     # Stage-A fast pre-filter residues
STAGE_A_CUTOFF    = 20.0   # Å — pairs above this skip full Kabsch

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H", "SEC": "U",
}

# ── Module-level globals shared with worker processes via fork ────────────────
_g_seqs   = []
_g_coords = []
_g_aln    = None
_g_sup    = None


def _worker_init():
    global _g_aln, _g_sup
    _g_aln = PairwiseAligner()
    _g_aln.mode = "global"
    _g_aln.substitution_matrix = substitution_matrices.load("BLOSUM62")
    _g_aln.open_gap_score   = -10
    _g_aln.extend_gap_score = -0.5
    _g_sup = SVDSuperimposer()


def _rmsd_worker(ij):
    """Two-stage RMSD: landmark pre-filter (Stage A) then full Kabsch (Stage B)."""
    i, j = ij
    aln = _g_aln.align(_g_seqs[i], _g_seqs[j])[0]
    idx_i, idx_j = aln.indices
    mask = (idx_i != -1) & (idx_j != -1)
    n_common = int(mask.sum())
    if n_common < 3:
        return 99.9

    ci_all = _g_coords[i][idx_i[mask]]
    cj_all = _g_coords[j][idx_j[mask]]

    # ── Stage A: Kabsch on N_LANDMARKS evenly-spaced residues ────────────────
    n_land   = min(N_LANDMARKS, n_common)
    land_idx = np.linspace(0, n_common - 1, n_land, dtype=int)
    _g_sup.set(ci_all[land_idx], cj_all[land_idx])
    _g_sup.run()
    if _g_sup.get_rms() > STAGE_A_CUTOFF:
        return 99.9                # too different — skip full alignment

    # ── Stage B: full Kabsch on all common residues ───────────────────────────
    _g_sup.set(ci_all, cj_all)
    _g_sup.run()
    return float(_g_sup.get_rms())


def get_sequence(struct):
    return "".join(_AA3TO1.get(r["res_name"], "X") for r in struct["ca_trace"])


_aligner_main = PairwiseAligner()
_aligner_main.mode = "global"
_aligner_main.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner_main.open_gap_score   = -10
_aligner_main.extend_gap_score = -0.5


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


if __name__ == "__main__":
    t_start = time.time()

    # ── 1. Load structures, skip short-chain flagged ──────────────────────────
    short_chain_ids = set()
    if os.path.exists(FLAG_FILE):
        with open(FLAG_FILE) as f:
            flag_data = json.load(f)
        short_chain_ids = set(flag_data.get("flagged", {}).keys())
        print(f"Short-chain flags: {len(short_chain_ids)} structures will be skipped")

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
    t_load = time.time()

    # ── 2. Align all structures to reference ──────────────────────────────────
    print(f"Aligning {N} structures to {ref_id}…")
    ref      = structures[ref_id]
    align_log = {ref_id: {"n_aligned": ref["residue_count"], "rmsd_to_ref": 0.0}}

    for i, pid in enumerate(pdb_ids):
        if pid == ref_id:
            continue
        n_aln, rmsd = kabsch_align(structures[pid], ref)
        align_log[pid] = {"n_aligned": n_aln, "rmsd_to_ref": round(rmsd, 3)}
        if (i + 1) % 100 == 0:
            print(f"  aligned {i + 1}/{N}")

    print(f"  aligned {N}/{N}  ✓")
    t_align = time.time()

    # ── 3. Flag RMSD outliers (> 3× median rmsd_to_ref) ──────────────────────
    rmsd_to_ref    = np.array([align_log[p]["rmsd_to_ref"] for p in pdb_ids if p != ref_id])
    median_rmsd    = float(np.median(rmsd_to_ref))
    outlier_cutoff = 3.0 * median_rmsd
    rmsd_outliers  = {
        pid for pid in pdb_ids
        if pid != ref_id and align_log[pid]["rmsd_to_ref"] > outlier_cutoff
    }
    print(f"\nRMSD-to-ref: median={median_rmsd:.2f} Å  outlier cutoff (3×)={outlier_cutoff:.2f} Å")
    print(f"  RMSD outliers excluded: {len(rmsd_outliers)}")
    if rmsd_outliers:
        for pid in sorted(rmsd_outliers):
            print(f"    {pid:<8} rmsd={align_log[pid]['rmsd_to_ref']:.2f} Å")

    active_ids = [p for p in pdb_ids if p not in rmsd_outliers]
    A          = len(active_ids)
    n_pairs    = A * (A - 1) // 2
    n_cpu      = os.cpu_count() or 1
    print(f"\nActive structures for clustering: {A}  ({n_pairs:,} pairs)  using {n_cpu} cores")

    # ── 4. Two-stage parallel pairwise RMSD ───────────────────────────────────
    print(f"Stage A cutoff: {STAGE_A_CUTOFF} Å on {N_LANDMARKS} landmarks")
    print("Computing pairwise RMSD matrix…")

    _g_seqs[:]   = [get_sequence(structures[p]) for p in active_ids]
    _g_coords[:] = [
        np.array([r["coords"] for r in structures[p]["ca_trace"]], dtype=np.float64)
        for p in active_ids
    ]

    pairs     = [(i, j) for i in range(A) for j in range(i + 1, A)]
    chunksize = max(1, len(pairs) // (n_cpu * 16))

    with Pool(processes=n_cpu, initializer=_worker_init) as pool:
        flat = pool.map(_rmsd_worker, pairs, chunksize=chunksize)

    t_rmsd = time.time()

    # count how many pairs were skipped by Stage A
    skipped_stageA = sum(1 for v in flat if v >= 99.9)
    print(f"  done — Stage A skipped {skipped_stageA:,}/{n_pairs:,} pairs "
          f"({100*skipped_stageA/max(n_pairs,1):.1f}%)")

    rmsd_matrix = np.zeros((A, A))
    for (i, j), v in zip(pairs, flat):
        rmsd_matrix[i, j] = rmsd_matrix[j, i] = v

    active_rmsd = rmsd_matrix[rmsd_matrix > 0]
    print(f"  min={active_rmsd.min():.2f}  max={active_rmsd.max():.2f}  "
          f"mean={active_rmsd[active_rmsd < 99.9].mean():.2f} Å")

    # ── 5. UPGMA clustering ───────────────────────────────────────────────────
    condensed = squareform(rmsd_matrix)
    Z         = linkage(condensed, method="average")
    labels    = fcluster(Z, t=CLUSTER_THRESHOLD, criterion="distance")

    clusters: dict[int, list[str]] = {}
    for pid, lbl in zip(active_ids, labels):
        clusters.setdefault(int(lbl), []).append(pid)

    print(f"\nClustering (UPGMA, {CLUSTER_THRESHOLD} Å):  {len(clusters)} clusters")
    for lbl in sorted(clusters):
        print(f"  C{lbl:<3}  {len(clusters[lbl]):>3} members:  "
              + ", ".join(clusters[lbl][:8]) + ("…" if len(clusters[lbl]) > 8 else ""))

    # ── 6. Save clusters.json ─────────────────────────────────────────────────
    condensed_list = [[round(v, 3) for v in row] for row in rmsd_matrix.tolist()]
    out = {
        "reference":                   ref_id,
        "cluster_threshold_angstroms": CLUSTER_THRESHOLD,
        "n_structures":                A,
        "n_clusters":                  len(clusters),
        "pdb_ids":                     active_ids,
        "rmsd_matrix":                 condensed_list,
        "cluster_assignments":         {pid: int(lbl) for pid, lbl in zip(active_ids, labels)},
        "clusters":                    {str(k): v for k, v in sorted(clusters.items())},
        "alignment_to_reference":      align_log,
        "excluded": {
            "short_chain":   sorted(short_chain_ids),
            "rmsd_outliers": sorted(rmsd_outliers),
        },
        "two_stage_rmsd": {
            "n_landmarks":    N_LANDMARKS,
            "stage_a_cutoff": STAGE_A_CUTOFF,
            "pairs_skipped":  skipped_stageA,
            "pairs_total":    n_pairs,
        },
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → {OUT_FILE}")

    t_end = time.time()
    print(f"\n{'─'*52}")
    print(f"  Load:      {t_load   - t_start:.1f}s")
    print(f"  Align:     {t_align  - t_load:.1f}s")
    print(f"  RMSD:      {t_rmsd   - t_align:.1f}s")
    print(f"  Cluster:   {t_end    - t_rmsd:.1f}s")
    print(f"  Total:     {t_end    - t_start:.1f}s")
