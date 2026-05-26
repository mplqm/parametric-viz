# viz3d Pipeline Skill

A reproducible pipeline for downloading a family of related structures from RCSB PDB,
extracting backbone geometry, clustering by structural similarity, generating smooth
interpolated transition paths between cluster centroids, and rendering an interactive
3-D viewer in the browser with no build step.

The reference run uses adenylate kinase (AK). Every design decision is documented so
the pipeline can be ported to a new protein family — or adapted to a different
scientific domain entirely — without guesswork.

---

## 1. Purpose

**Use this pipeline when you want to:**

- Compare a family of related 3-D structures (proteins, RNA, small molecules, …)
- Identify structural clusters within that family
- Visualise the conformational "landscape" — what does the space between clusters
  look like as a smooth animation?
- Deliver a self-contained, shareable browser viewer that needs only
  `python3 -m http.server` to run

**This is not the right tool when you need:**

- Full-atom molecular dynamics (use GROMACS / AMBER)
- Crystallographic refinement or electron-density analysis
- Sub-Å precision structural comparison (use dedicated tools like TM-align, FATCAT)

---

## 2. Prerequisites

### Python packages

```
python3-biopython   # BioPython 1.80+ (apt: sudo apt install python3-biopython)
python3-numpy       # installed automatically as biopython dependency
python3-scipy       # sudo apt install python3-scipy
```

No virtualenv or conda needed on Debian/Ubuntu — system packages are sufficient.

### Folder layout expected at the start

```
project/
├── hello_pdb.py            ← proof-of-concept / sanity check (optional)
├── download_family.py
├── extract_geometry.py
├── align_and_cluster.py
├── build_transitions.py
└── index.html
```

The pipeline creates these directories automatically:

```
project/
├── pdb_files/              ← raw PDB downloads  (~12 MB for 20 AK structures)
├── geometry/               ← per-structure CA-trace JSON  (~1 MB)
├── clusters.json           ← RMSD matrix + cluster assignments  (~8 KB)
└── transitions.json        ← interpolated frames for all cluster pairs  (~12 MB)
```

### Hardware

Any modern laptop. The slowest step (RMSD matrix computation, 190 pairwise
BLOSUM62 alignments + SVD) takes ~2 minutes on a single CPU core for 20 structures.
The browser viewer runs comfortably on integrated graphics at 300 Cα atoms per frame.

---

## 3. Step-by-Step Pipeline

### Script 1 — `hello_pdb.py`  *(sanity check, optional)*

**What it does:** Downloads one PDB file and counts atoms.  
**Run:** `python3 hello_pdb.py`  
**Inputs:** Hardcoded PDB ID (`4AKE`)  
**Outputs:** `4AKE.pdb` in the working directory, atom count printed  
**Adjust:** Change `PDB_ID` to any 4-character RCSB accession  

---

### Script 2 — `download_family.py`

**What it does:**
1. Queries the RCSB Search API v2 for structures matching a search term
2. Filters to X-ray crystallography structures only (which have a resolution value)
3. Sorts by `resolution_combined` ascending (best resolution first)
4. Downloads the top N PDB files into `pdb_files/`
5. Prints PDB ID, resolution (Å), and atom count for each

**Run:** `python3 download_family.py`  
**Inputs:** RCSB Search API (internet required)  
**Outputs:** `pdb_files/*.pdb` (20 files, ~600 KB each)  

**Key parameters to adjust:**

| Parameter | Location | Default | Effect |
|-----------|----------|---------|--------|
| Search term | `"value": "adenylate kinase"` | `"adenylate kinase"` | Change to your protein family |
| Experimental method | `"value": "X-ray"` | `"X-ray"` | `"NMR"` has no resolution; `"EM"` for cryo-EM |
| Number of structures | `"rows": 20` | `20` | Increase for broader sampling; >50 makes the RMSD matrix slow |
| Sort field | `sort_by: rcsb_entry_info.resolution_combined` | resolution | Could sort by deposition date, R-factor, etc. |

**RCSB API note:** The experimental method filter value must be `"X-ray"` (not
`"X-RAY DIFFRACTION"` — case and exact string matter).

---

### Script 3 — `extract_geometry.py`

**What it does:**
1. Reads every `.pdb` file in `pdb_files/`
2. Extracts Cα (alpha-carbon) atoms only — one per residue, first model only
3. Skips HETATM records (ligands, water) and non-standard residues lacking Cα
4. Computes the centroid of the Cα cloud and subtracts it (centres at origin)
5. Writes one JSON per structure to `geometry/`

**Run:** `python3 extract_geometry.py`  
**Inputs:** `pdb_files/*.pdb`  
**Outputs:** `geometry/<PDB_ID>.json`  

**Output fields per file:**
```json
{
  "pdb_id": "3QMA",
  "chains": ["A"],
  "centroid_original": [-6.7, -0.9, 12.7],
  "residue_count": 145,
  "ca_trace": [
    { "chain": "A", "res_seq": 1, "res_name": "MET", "coords": [-1.8, -16.0, 3.3] },
    ...
  ]
}
```

**Key parameters to adjust:**

- **`PDB_DIR`** / **`OUT_DIR`** — change folder names if needed
- **For non-protein structures** — replace the Cα selection with your atom type
  (see Domain Adaptation section)

---

### Script 4 — `align_and_cluster.py`

**What it does:**
1. Loads all geometry JSONs
2. For each structure, runs a BLOSUM62 global pairwise sequence alignment against
   the reference structure to establish residue correspondence
3. Applies Kabsch (SVD least-squares) superimposition of each structure onto the
   reference, updating coordinates in memory
4. Computes an N×N pairwise RMSD matrix using Kabsch alignment on
   sequence-matched Cα pairs
5. Runs UPGMA agglomerative hierarchical clustering on the RMSD matrix
6. Cuts the dendrogram at `CLUSTER_THRESHOLD` Å
7. Writes `clusters.json`

**Run:** `python3 align_and_cluster.py`  
**Inputs:** `geometry/*.json`  
**Outputs:** `clusters.json`  

**Key parameters to adjust:**

| Parameter | Default | When to change |
|-----------|---------|----------------|
| `CLUSTER_THRESHOLD` | `15.0` Å | Lower (e.g. 3 Å) for a highly conserved family; higher (e.g. 20 Å) for a very diverse family. Start at 3 Å, raise until you get 3–8 biologically meaningful clusters |
| `method="average"` in `linkage()` | UPGMA | `"ward"` for compact spherical clusters; `"complete"` for conservative merging |
| BLOSUM62 gap penalties | open=−10, extend=−0.5 | Widen gaps (open=−6) for distantly related sequences; tighten for closely related |

**Choosing the reference structure:**
The script uses the alphabetically first structure if `4AKE` is not present.
For your own dataset, explicitly set `ref_id` to a well-resolved, full-length,
representative structure — not a truncated construct or an outlier.

---

### Script 5 — `build_transitions.py`

**What it does:**
1. Re-runs the reference alignment (same code as Script 4) to get all structures
   into the reference coordinate frame
2. For each cluster, maps every member's Cα atoms to reference sequence positions
   (via BLOSUM62 alignment) and averages coordinates → cluster centroid
3. For every pair of cluster centroids, finds the residues present in both
4. Generates `N_FRAMES` interpolated frames between the two centroids using
   smoothstep easing
5. Writes `transitions.json`

**Run:** `python3 build_transitions.py`  
**Inputs:** `geometry/*.json`, `clusters.json`  
**Outputs:** `transitions.json`  

**Key parameters to adjust:**

| Parameter | Default | Effect |
|-----------|---------|--------|
| `N_FRAMES` | `60` | Frames per transition. 30 = snappier; 120 = smoother slow-motion |
| `CLUSTER_THRESHOLD` (in align_and_cluster.py) | `15.0` Å | Fewer clusters → fewer transitions; more clusters → more transitions but smaller JSON |

**Note on centroid residue counts:** A multi-member cluster centroid only includes
residue positions present in **all** members. Diverse clusters (different species,
truncated constructs) will have fewer centroid residues. Singleton clusters (one
structure) keep all their residues.

---

### Script 6 — `index.html`

**What it does:** Self-contained Three.js viewer — no npm, no build step.  
**Serve:** `python3 -m http.server 8080` then open `http://localhost:8080`  
**Inputs:** `transitions.json` (fetched at runtime)  

Controls:
- **Dropdown** — select any of the pairwise transitions
- **Slider** — scrub frames 0–59 manually; dragging pauses auto-play
- **▶ Play / ⏸ Pause** — auto-advance at selected speed
- **Speed selector** — ½×, 1×, 2×, 4× (6 / 12 / 24 / 48 fps)
- **↩ Reset** — return to frame 0
- **Mouse** — left-drag rotates, right-drag pans, scroll zooms (OrbitControls)

---

## 4. Domain Adaptation Guide

The pipeline has three swappable layers:

```
[Data source]  →  [Geometry extraction]  →  [Alignment + clustering]  →  [Viewer]
```

Only the first two layers are domain-specific. The clustering and viewer are
generic once you have `geometry/*.json` files with a `ca_trace`-style list
of 3-D points.

### Adapting to a different protein family

Change only the search term in `download_family.py`:

```python
"value": "lysozyme"          # or "hsp90", "GPCR", "beta-lactamase", …
```

Everything else works unchanged.

### Adapting to RNA structures

RNA uses **C3′** (3-prime carbon) as the backbone trace atom, equivalent to
Cα in proteins. In `extract_geometry.py`, change the atom selection:

```python
# protein
if "CA" not in residue: continue
ca = residue["CA"]

# RNA
if "C3'" not in residue: continue
ca = residue["C3'"]
```

For alignment in `align_and_cluster.py`, replace BLOSUM62 with a nucleotide
substitution matrix (or use simple match/mismatch scoring since there are only
4 nucleotides):

```python
_aligner.match_score    =  2
_aligner.mismatch_score = -1
```

### Adapting to small molecules / ligands

There is no backbone trace atom. Options:

1. **All heavy atoms:** Use every non-hydrogen atom. Skip the sequence alignment
   entirely; match atoms by atom name (`C1`, `N2`, etc.) or by index if the
   molecules are identical.
2. **Pharmacophore points:** Extract centroid of each ring system, each H-bond
   donor/acceptor. This reduces a 50-atom molecule to 5–10 meaningful points.

In `align_and_cluster.py`, replace the BLOSUM62 sequence-matching block with
direct positional matching (same atom names in both structures).

### Adapting to cryo-EM structures

CryoEM structures in RCSB are deposited as PDB/mmCIF files with atomic models —
they work identically to X-ray structures. Change the experimental method filter:

```python
"value": "ELECTRON MICROSCOPY"
```

Resolution for cryo-EM is stored in the same field (`rcsb_entry_info.resolution_combined`),
so the sort and download logic is unchanged.

### Adapting to non-biological 3-D data (e.g. materials, astronomy)

Skip Scripts 1–3 entirely. Write your own script that produces the same
`geometry/<ID>.json` format (see Output Schema section). As long as each file
has a list of 3-D points in `ca_trace[].coords`, Scripts 4–6 will work without
modification.

For alignment, replace BLOSUM62 with identity matching or pure ICP
(Iterative Closest Point) if there is no sequence concept.

---

## 5. Key Decisions Documented

### Why BLOSUM62 for residue matching?

Residue numbers in PDB files are assigned by the depositor and are inconsistent
across structures from different organisms or labs. Simple number matching fails:
`2BZZ` numbers its residues `1000–1134`. Chain labels also vary (`X` instead of `A`).

BLOSUM62 global sequence alignment solves this by finding the best residue
correspondence based on evolutionary substitution probabilities, regardless of
numbering scheme. It handles insertions and deletions cleanly and is robust for
sequences with 30–100% identity — the typical range within a protein family.

Alternatives considered:
- **Residue number matching:** Fails for non-standard numbering (as seen with 2BZZ)
- **Positional matching (first N residues):** Fails when structures have different
  N-terminal truncations
- **Smith-Waterman local alignment:** Better for very distantly related sequences,
  but overkill for a family search that already filters by name

### Why 15 Å clustering threshold?

At 3 Å (strict), the AK dataset yields 18–19 singleton clusters because the
structures span genuinely large conformational changes (7–10 Å RMSD for the LID
domain alone between open and closed states) and evolutionary diversity across
organisms. Nearly every structure would be its own cluster — not useful.

At 15 Å, UPGMA agglomeration produces 13 clusters including 5 multi-member groups
that reflect real biological groupings (near-identical depositions, same organism
different conditions, short constructs from one species).

**General rule for choosing the threshold:**
1. Run `align_and_cluster.py` at `CLUSTER_THRESHOLD = 3.0`
2. Count singleton clusters. If >60% of structures are singletons, raise the threshold
3. Aim for 3–8 clusters with at least 2–3 multi-member clusters
4. Inspect whether the multi-member clusters make biological sense

### Why Cα atoms only?

- **Backbone trace:** Cα positions capture overall fold. Side-chain positions vary
  with crystal contacts and local flexibility; they add noise to clustering.
- **Performance:** ~1/10 the atom count of the full structure. A 300-residue protein
  has ~300 Cα vs ~2400 heavy atoms.
- **Universality:** Every standard amino acid has exactly one Cα, so the trace is
  unambiguous and one-dimensional (residue index → 3-D point).
- **Alignment:** Cα-only RMSD is the standard metric in structural biology. Papers
  report "Cα RMSD" by convention; our values are directly comparable.

### Why smoothstep easing for interpolation?

Linear interpolation gives constant velocity: the transition starts and ends
abruptly, which looks mechanical in an animation.

Smoothstep: `f(t) = t²(3 − 2t)` has zero first derivative at `t=0` and `t=1`,
meaning the motion accelerates from rest and decelerates to rest. This feels
more physical and is easier to follow visually.

For higher quality, cubic or Catmull-Rom spline interpolation through all cluster
centroids would give a continuous path. Smoothstep per pair is sufficient for an
exploratory viewer and requires no global path planning.

### Why HTTP server instead of `file://`?

Modern browsers enforce the **same-origin policy** for `fetch()`. A page loaded
via `file://` cannot fetch another local file — the browser treats them as
cross-origin and blocks the request. `python3 -m http.server` serves all files
from the same origin (`http://localhost:8080`) so `fetch('transitions.json')`
succeeds.

Alternative: embed `transitions.json` as a JavaScript variable inside `index.html`
(`const DATA = { ... };`). This eliminates the server requirement but makes
`index.html` very large (~12 MB) and harder to update independently.

---

## 6. Known Edge Cases

### Tetramers and large multi-chain assemblies (e.g. `3U7Q`)

`3U7Q` is a tetrameric adenylate kinase (chains A–D, 1998 Cα atoms). The pipeline
includes all chains from the asymmetric unit. Consequences:

- The centroid of a multi-chain cluster will have far more residues than a monomer,
  and pairwise RMSD to monomers will be very high (28–40 Å in the AK run)
- `3U7Q` correctly ends up as a singleton cluster
- In the viewer it appears as 4 superimposed chains, which can look like a dense blob

**Fix:** In `extract_geometry.py`, add a chain filter to keep only chain A:
```python
for chain in model:
    if chain.id != 'A': continue   # add this line
```
Or filter at download time by choosing structures with `assembly_count == 1`
using the RCSB data API.

### Non-standard residue numbering (e.g. `2BZZ`: residues 1000–1134)

Some depositors use non-contiguous or offset numbering (domain constructs,
engineered variants, numbering relative to a larger parent protein).
Residue-number matching silently fails for these — the fallback to positional
matching pairs the wrong residues and gives high RMSD.

The BLOSUM62 sequence alignment in `align_and_cluster.py` handles this correctly
because it matches by amino acid identity, not residue number. If you use a
simpler matching strategy, always validate with a spot-check:
`python3 -c "import json; d=json.load(open('geometry/2BZZ.json')); print(d['ca_trace'][0])"`.

### Non-standard chain labels (e.g. `2AT3`: chain `X`)

PDB chain labels are depositor-assigned and can be any single character. Matching
by `(chain, res_seq)` fails when the same protein is deposited as chain `A` in
one structure and chain `X` in another. The BLOSUM62 alignment in Script 4 handles
this; simple res_seq-only matching partially handles it; positional matching
ignores it.

### Truncated constructs (e.g. `6MM2`: 99 residues, `4AFF`: 110 residues)

Short constructs cover only part of the reference protein. They will:
- Cluster together if they cover the same domain (4AFF and 6MM2 cluster at ~5 Å)
- Produce centroids with few residues (centroid of cluster 9 has only 13 Cα)
- Generate transitions with few tracked residues (5–13 residues for cluster 9 pairs)

Transitions with fewer than ~20 residues are structurally uninformative. Consider
filtering them out of the viewer dropdown, or applying a minimum-residue threshold:
```python
if t['n_residues'] < 20:
    continue  # skip in build_transitions.py
```

### Near-identical duplicate depositions (e.g. `5ONK` / `6HRI`: RMSD 0.04 Å)

Some structures are re-deposited refinements of the same crystal. They produce a
2-member cluster with near-zero RMSD and a trivial transition (barely any motion).
These are valid but visually uninteresting. They are correctly grouped and can
be left in the dataset.

---

## 7. Output Schema

### `geometry/<PDB_ID>.json`

```json
{
  "pdb_id":             "3QMA",
  "chains":             ["A"],
  "centroid_original":  [-6.69, -0.95, 12.75],
  "residue_count":      145,
  "ca_trace": [
    {
      "chain":    "A",
      "res_seq":  1,
      "res_name": "MET",
      "coords":   [-1.81, -16.03,  3.30]
    }
  ]
}
```

- `coords` are in Ångströms, centred at the origin (centroid subtracted)
- One entry per residue; HETATM and water excluded

### `clusters.json`

```json
{
  "reference":                   "1OD8",
  "cluster_threshold_angstroms": 15.0,
  "pdb_ids":                     ["1OD8", "2AT3", ...],
  "rmsd_matrix":                 [[0.0, 21.2, ...], ...],
  "n_clusters":                  13,
  "cluster_assignments":         {"1OD8": 3, "2AT3": 8, ...},
  "clusters":                    {"1": ["6KFN"], "5": ["2BZZ","5ONK","6HRI"], ...},
  "alignment_to_reference":      {
    "3QMA": { "n_aligned": 143, "rmsd_to_ref": 19.193 }
  }
}
```

- `rmsd_matrix[i][j]` is the Kabsch RMSD between structures `pdb_ids[i]` and
  `pdb_ids[j]` in Ångströms, on BLOSUM62-matched Cα atoms
- `cluster_assignments` maps each PDB ID to its integer cluster label

### `transitions.json`

```json
{
  "reference":     "1OD8",
  "n_clusters":    13,
  "n_frames":      60,
  "interpolation": "smoothstep",
  "centroids": {
    "1": {
      "members":    ["6KFN"],
      "n_residues": 270,
      "res_seq":    [2, 3, 4, ...],
      "coords":     [[-24.7, -8.1, 8.6], ...]
    }
  },
  "transitions": [
    {
      "from_cluster": "1",
      "to_cluster":   "2",
      "from_members": ["6KFN"],
      "to_members":   ["7NIY"],
      "n_residues":   180,
      "res_seq":      [2, 3, 4, ...],
      "frames":       [
        [[-24.7, -8.1, 8.6], ...],
        ...
      ]
    }
  ]
}
```

- `centroids[label].coords[i]` is the average Cα position (in reference frame,
  Å) for the residue `res_seq[i]`, averaged over all cluster members
- `transitions[t].frames[k][i]` is the 3-D position `[x, y, z]` of residue
  `res_seq[i]` at frame `k` (0 = from-centroid, 59 = to-centroid)
- Frames use smoothstep interpolation: `f(k/59) = t²(3−2t)` applied to each axis
- The file uses compact JSON (`separators=(",",":")`) — no whitespace — to keep
  size manageable (~12 MB for 76 transitions × 60 frames × ~150 residues)

---

## 8. Viewer Notes

### Three.js setup

The viewer uses Three.js 0.165.0 loaded via an **importmap** — no npm, no bundler:

```html
<script type="importmap">
{
  "imports": {
    "three":         "https://cdn.jsdelivr.net/npm/three@0.165.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.165.0/examples/jsm/"
  }
}
</script>
<script type="module">
  import * as THREE from 'three';
  import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
  import { Line2 } from 'three/addons/lines/Line2.js';
```

Importmaps are supported in Chrome 89+, Firefox 108+, Safari 16.4+. For older
browsers, replace with a bundled build or use three.js r128 with the old UMD format.

### Rendering architecture

| Object | Three.js type | Why |
|--------|---------------|-----|
| Cα spheres | `InstancedMesh` | Single draw call for up to 2000 spheres; positions updated in-place each frame via `setMatrixAt()` |
| Backbone line | `Line2` (from addons) | Supports pixel-width lines (unlike `LineBasicMaterial` which is always 1 px in WebGL); vertex colours for N→C gradient |

### Frame update cost

Each frame change does:
1. N calls to `setMatrixAt()` + one `instanceMatrix.needsUpdate = true`
2. One `LineGeometry.setPositions()` call per backbone segment

For 300 residues at 48 fps this is well within browser budget (<2 ms per frame).
Performance only becomes a concern above ~5000 residues per frame.

### Colour scheme

Residues are coloured by their `res_seq` value normalised to the global min/max
across all transitions, mapped to HSL hue 0°→270° (red N-terminal → violet
C-terminal). This is consistent across all transitions so the same residue always
has the same colour regardless of which transition is shown.

### Backbone gap detection

The backbone line is broken into segments wherever consecutive `res_seq` values
differ by more than 6. This prevents false connections across chain breaks or
missing residue ranges. Each segment becomes an independent `Line2` object.

### Serving the viewer

```bash
python3 -m http.server 8080
# open http://localhost:8080
```

For production sharing:
- GitHub Pages works directly (push the repo, enable Pages on `main`)
- The 12 MB `transitions.json` may load slowly on slow connections;
  split by transition group or compress with gzip if needed
  (`Content-Encoding: gzip` is supported by most static hosts)

---

## 9. Worked Example — Adenylate Kinase Reference Run

### 1. Download

```bash
python3 download_family.py
```

Expected: 20 PDB files in `pdb_files/`. Best resolution: **3IP0** and **6KFN**
at 0.89 Å. Note that `4AKE` (the canonical reference, 2.0 Å) does **not** appear
in the top-20-by-resolution set.

### 2. Extract geometry

```bash
python3 extract_geometry.py
```

Expected output:
```
PDB ID   Chains   Residues
--------------------------
1OD8     A             301
2AT3     X             184
2BZZ     A             135     ← residues 1000–1134 (non-standard numbering)
...
3U7Q     A,B,C,D      1998     ← tetramer
...
6MM2     A              99     ← truncated construct
```

### 3. Align and cluster (15 Å threshold)

```bash
python3 align_and_cluster.py
```

Expected: 13 clusters including:
- **C5** `{2BZZ, 5ONK, 6HRI}` — 5ONK/6HRI are near-identical (0.04 Å); 2BZZ joins at ~13 Å
- **C7** `{3QMA, 8ESU, 9V14}` — 3QMA/8ESU very similar (2.06 Å)
- **C9** `{3IP0, 4AFF, 6MM2}` — three short/truncated constructs
- **C13** `{3U7Q}` — tetramer outlier, RMSD > 28 Å to everything else

RMSD matrix highlights:
- Most monomer pairs: **13–25 Å** (genuine structural diversity; AK has large
  domain motions and spans bacteria→eukaryotes)
- 3QMA/8ESU: **2.06 Å** (same protein, same state, slightly different conditions)
- 5ONK/6HRI: **0.04 Å** (essentially identical depositions)
- 3U7Q to any monomer: **28–40 Å** (tetramer vs monomer is not a fair comparison)

### 4. Build transitions

```bash
python3 build_transitions.py
```

Expected:
- **78 potential pairs**, 2 skipped (C5↔C9 and C7↔C9 share no common residues
  after BLOSUM62 alignment — the short constructs in C9 and the partial coverage
  of C5's centroid do not overlap)
- **76 transitions** written
- `transitions.json`: ~12 MB
- Centroid residue counts range from 13 (C9, 3 truncated constructs) to 301 (C3, 1OD8 alone)

### 5. View

```bash
python3 -m http.server 8080
```

Open `http://localhost:8080`. Visually interesting transitions to explore:

| Transition | What you see |
|------------|-------------|
| C3→C5 (1OD8 → {2BZZ,5ONK,6HRI}) | Full-length reference morphs to a mid-protein fragment |
| C7→C8 ({3QMA,8ESU,9V14} → {2AT3,7BCU}) | Similar-length structures, different organisms |
| C3→C13 (1OD8 → 3U7Q tetramer) | Dramatic: monomer centroid morphs toward a tetramer sub-unit position |

The colour gradient (red N-term → violet C-term) makes it easy to follow which
domain is moving during a transition. The LID domain (C-terminal, violet end)
shows the largest displacements in most AK transitions, consistent with the
known open/closed mechanism.
