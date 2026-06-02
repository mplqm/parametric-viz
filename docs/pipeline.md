# Pipeline Documentation

How the parametric-viz pipeline works, and how to adapt it for a new domain.

Built by one person and Claude. Open source under MIT.

---

## What the pipeline does

Takes publicly available scientific structure data and produces:

- A set of geometry JSONs — one per structure, in a common reference frame
- A cluster assignment file — which structures are conformationally similar
- A set of transition JSONs — smooth interpolations between cluster centroids
- A viewer — WebXR AR on Android, AR Quick Look on iOS, Three.js on desktop

The pipeline is domain-agnostic. It was built on proteins but the architecture assumes nothing about the underlying science. Swap the download script and geometry extractor, and it runs on crystals, neurons, astronomical objects, or anything else with 3D structure data in a public database.

---

## Architecture overview

```
Public database
      ↓
download_family.py        — fetch raw structure files
      ↓
extract_geometry.py       — parse → 3D coordinates → centre → JSON
      ↓
align_and_cluster.py      — pairwise RMSD → cluster → keyframe states
      ↓
build_transitions.py      — interpolate between states → transition JSONs
      ↓
index.html                — load on demand → render in AR or desktop
      ↓
publish.py                — Zenodo upload → DOI → patch JSONs → GitHub issues → repo public
```

Each stage is a standalone script. Run them in sequence. Outputs of each stage are inputs to the next. Nothing is stateful between runs — restart any stage without affecting others.

---

## Dependencies

```bash
pip install biopython numpy scipy requests tqdm --break-system-packages
```

- **BioPython** — SVDSuperimposer for Kabsch/SVD alignment, PDB parsing
- **NumPy** — coordinate arithmetic, matrix operations
- **SciPy** — agglomerative clustering
- **requests** — database API calls
- **tqdm** — progress bars for long runs

Python 3.9+ required. Tested on 3.11.2.

---

## Stage 1 — Download (`download_family.py`)

Queries a public database API, filters by quality, downloads raw structure files.

**Key decisions:**
- Filter threshold is set per domain — proteins use resolution ≤ 2.5Å
- Skip structures already on disk — safe to rerun after interruption
- Never relax quality filters to hit round numbers
- Log all exclusions with reasons

**Output:** `pdb_files/` — raw files, gitignored. Always re-downloadable. Never commit these.

**To adapt for a new domain:**
Replace the API query and download logic. Everything else stays the same. The output just needs to be a folder of raw files that `extract_geometry.py` can read.

---

## Stage 2 — Geometry Extraction (`extract_geometry.py`)

Parses raw structure files, extracts 3D coordinates, centres them, saves one JSON per structure.

**What gets extracted:**

| Domain | Coordinates extracted |
|---|---|
| Proteins | CA backbone atoms — one per residue |
| Crystals | Unit cell atomic positions |
| Neurons | Dendritic trace coordinates |
| Astronomy | Stellar parameter positions in parameter space |

**Centring — the most important invariant:**

Every structure is centred at the combined centroid of ALL chains together. Never per-chain. This places all structures in the same reference frame so RMSD comparisons are valid.

```python
# CORRECT
all_coords = np.vstack([coords_per_chain[c] for c in chains])
centroid = all_coords.mean(axis=0)
centered_coords = raw_coords - centroid

# WRONG — breaks the shared coordinate space
centroid = chain_coords.mean(axis=0)  # per-chain centring
```

This was verified against 10JX.json (2 chains, 282 residues). Do not change it.

**Quality flags:**
- Chain A < 50% of median chain length → skip + log
- Centring failure → fix before proceeding

**Output:** `geometry/` — one JSON per structure. Committed to GitHub.

**To adapt for a new domain:**
Change the coordinate extraction function. The centring logic, JSON output format, and quality flags stay identical.

---

## Stage 3 — RMSD Matrix + Clustering (`align_and_cluster.py`)

Computes pairwise structural similarity across all structures, then groups them into conformational clusters.

### RMSD — two-stage approach

Computing full pairwise RMSD across thousands of structures is O(n²). The two-stage approach makes this tractable:

**Stage A — Landmark pre-filter:**
Select 20 evenly-spaced residues as landmarks. Compute rough RMSD on landmarks only. If landmark RMSD exceeds a threshold, skip the full computation. Eliminates ~80% of pairs.

**Stage B — Full Kabsch/SVD:**
For remaining pairs, run full alignment via BioPython's SVDSuperimposer. This is the correct algorithm — it finds the optimal rotation to minimise RMSD between two structures.

```python
from Bio.SVDSuperimposer import SVDSuperimposer

sup = SVDSuperimposer()
sup.set(reference_coords, mobile_coords)
sup.run()
rmsd = sup.get_rms()
```

**Residue matching:**
Structures from the same protein family often have non-standard residue numbering. BLOSUM62 sequence alignment is used to identify which residues correspond across structures before computing RMSD.

**Parallelisation:**
The RMSD matrix is computed in parallel across all available CPU cores:

```python
from multiprocessing import Pool

# Always guard with this on Windows
if __name__ == '__main__':
    with Pool(processes=8) as pool:
        results = pool.map(compute_rmsd_pair, pairs)
```

Scale `processes` to your machine's thread count.

**Outlier flagging:**
Structures with RMSD > 3× the median of the full set are flagged in metadata but retained in provenance. They're not deleted — a domain expert may find them interesting.

### Clustering

Agglomerative clustering at a 15Å threshold groups structures into conformational families. Each cluster centroid becomes a keyframe state in the viewer.

The 15Å threshold was chosen for proteins — it produces meaningful conformational groups without over-splitting. Adjust per domain if needed, but document any change.

**Output:** `clusters.json` — RMSD matrix, cluster assignments, centroid coordinates. Committed to GitHub.

---

## Stage 4 — Transition Generation (`build_transitions.py`)

Interpolates smoothly between every pair of cluster centroids. Each transition is a sequence of 60 frames that the viewer loads on demand.

**Interpolation:**
Linear interpolation between centroid coordinate sets, with smoothstep easing applied to the parameter t:

```python
def smoothstep(t):
    return t * t * (3 - 2 * t)

frames = []
for i in range(60):
    t = smoothstep(i / 59)
    frame = start_coords + t * (end_coords - start_coords)
    frames.append(frame.tolist())
```

Smoothstep gives zero velocity at start and end — transitions feel natural rather than mechanical.

**Skipping:**
Transitions between cluster pairs with no common residues are skipped entirely. Log the skip with both cluster IDs.

**Output:** `transitions/` — one JSON per transition. Gitignored. Zenodo only.

Transition files are large in aggregate (6GB+ for adenylate kinase alone). GitHub can't hold them. Zenodo gives them a DOI and permanent storage.

**Lazy loading:**
The viewer fetches transition files on demand — only when a user selects that transition. Nothing is preloaded. This keeps the viewer fast regardless of how many transitions exist.

---

## Viewer (`index.html`)

Single file. Three rendering tiers, one URL, auto-detected on load:

**Android Chrome** → WebXR AR
- Hit testing places protein on real surface
- Reticle shown before placement, confirmed on tap
- Pinch to scale, two-finger rotate
- Floating UI panel in AR space
- WebXR light estimation, fallback to ambient + directional
- Subsampled to 214 CA atoms max for mobile 60fps

**iOS Safari** → AR Quick Look via USDZExporter
- USDZ generated client-side from centroid geometry
- Triggered via `window.location.href` blob URL
- Scale baked into vertex positions (mesh.scale ignored by USDZExporter)
- Scale factor 0.0015 — maps 72Å → ~10cm
- Rainbow gradient via 256×1 canvas texture (vertex colours unsupported in AR Quick Look)
- Backbone as CylinderGeometry (LineSegments silently skipped)
- InstancedMesh converted to merged BufferGeometry before export
- Static conformation only — no animated transitions in ARKit

**Desktop** → Three.js with OrbitControls
- Fallback only
- Transition dropdown, frame slider, play/pause, speed control
- No DOF, no fog, no post-processing effects

---

## Forking for a New Domain

Minimum changes required:

1. **`download_family.py`** — replace API endpoint, query parameters, and file download logic. Output: a folder of raw structure files.

2. **`extract_geometry.py`** — replace the coordinate extraction function for your data format. Keep: centring logic, JSON output schema, quality flags.

3. **Quality threshold** — set the equivalent of "resolution ≤ 2.5Å" for your domain. Document it. Never change it mid-run.

4. **Clustering threshold** — 15Å works for proteins. Evaluate whether it produces meaningful groups for your domain. Document any change.

Everything else — RMSD computation, clustering, transition generation, viewer, JSON schema, git workflow, Zenodo publishing — runs unchanged.

---

## JSON Schema

Every domain produces identical JSON structure. See `docs/schema.md` for full specification.

Critical fields:
- `pipeline_version` — tracks which run generated each file
- `stage` — 1=automated, 2=expert validated
- `flagged` — triggers GitHub issue creation on publish
- `stage2_annotations` — null until a domain expert contributes

---

## Git Workflow

```bash
git add -A
git commit -m "Domain: [name] — [N] structures, [N] clusters, [N] transitions"
git push origin main
```

Credentials stored via `git config credential.helper store`. Never embed PAT in commands.

What goes to GitHub: pipeline scripts, geometry JSONs, clusters.json, README, SKILL.md, pipeline.md.

**Repo structure** (key files):
```
viz3d_pipeline/
  download_family.py
  extract_geometry.py
  align_and_cluster.py
  build_transitions.py
  publish.py                 — Zenodo + GitHub publish automation
  index.html
  domain_config/             — one JSON per domain, Zenodo metadata
    proteins.json
  geometry/
  clusters.json
  transitions/               — Zenodo only
  docs/
    pipeline.md
    schema.md
    contributing.md
    publish.md
```

What goes to Zenodo only: transitions/ (too large for GitHub, needs DOI).

What goes nowhere: pdb_files/ (always re-downloadable from source).

---

## Compute Notes

**Bottleneck:** RMSD matrix. Scales as O(n²) before the landmark pre-filter, O(0.2n²) after.

**Z420 baseline** (Xeon E5-1620 v2, 8 threads):
- 942 structures → ~443,000 pairs → ~8–12 hours with parallelisation

**At domain scale** (100,000+ structures), consider:
- Sub-family RMSD — run matrix within sub-families, not globally
- Cloud burst for the RMSD stage — a 96-core instance reduces weeks to hours
- The rest of the pipeline (extraction, transitions, viewer) scales linearly and runs fine on the Z420

---

## Adenylate Kinase — Reference Results

| Metric | Value |
|---|---|
| Structures processed | 942 |
| Quality filter pass rate | 94.2% |
| RMSD pairs computed | ~443,000 |
| Conformational clusters | 313 |
| Animated transitions | 30,735 |
| Transition data | 6.0 GB |
| Centring | Confirmed correct ✅ |
| iOS AR | Confirmed working May 2026 ✅ |

Use these as a sanity check when running a new domain. If cluster counts or transition volumes are wildly different in proportion, investigate before proceeding.

---

## Stage 5 — Publish (`publish.py`)

Automates the full Zenodo + GitHub publish sequence. Run after git commit. No manual steps required.

```bash
python3 publish.py --domain proteins --sandbox   # test run first
python3 publish.py --domain proteins             # production
```

**What it does in sequence:**

1. Uploads `transitions/` folder to Zenodo as a dataset release
2. Uploads `geometry/` folder to the same deposit
3. Sets metadata — domain, family, pipeline version, structure count
4. Publishes the deposit and retrieves the DOI
5. Patches `zenodo_doi` into all model JSONs
6. Commits the DOI-updated JSONs to git and pushes
7. Opens GitHub issues for all flagged structures (Stage 2 invitation)
8. Makes the GitHub repo public via GitHub API

**Requires two tokens** set as environment variables:

```bash
export ZENODO_TOKEN=your_zenodo_token
export GITHUB_TOKEN=your_github_pat
```

**Always run sandbox first.** Zenodo sandbox is a full dry-run environment — identical behaviour, no public deposit. Confirm the DOI comes back and JSONs update correctly before running production.

**Domain config** lives in `domain_config/<domain>.json` — title, description, keywords, creator metadata for the Zenodo record.

See `docs/publish.md` for full usage, token setup, and domain config format.

---

## Animation Hold (Pending ⬜)

Not yet implemented. Specified behaviour for desktop and WebXR AR viewers:

- 2 second pause at each end state before the next transition plays
- Gives the viewer time to read the structure before motion resumes
- Applies to both desktop Three.js and WebXR AR
- Does not apply to iOS AR Quick Look — ARKit shows static conformation only

Implementation: track animation state in the viewer loop; when a transition completes, set a hold timer before queuing the next transition.

---

*Built by one person and Claude, on a workstation from 2013.*
