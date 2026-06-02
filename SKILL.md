---
name: parametric-viz-pipeline
description: Run the full parametric-viz pipeline for a new scientific domain — download, geometry extraction, RMSD matrix, clustering, transition generation, git commit. Use this skill whenever the user asks to run the pipeline, process a new domain, add structures, generate transitions, or anything related to the parametric-viz scientific visualization pipeline. Also use when the user mentions domain names from the roadmap (astronomy, crystals, neurons, archaeology, climate) in the context of processing or visualization.
---

# Parametric Viz Pipeline

A fully automated pipeline that transforms publicly available scientific structure data into interactive parametric 3D visualizations. Domain-agnostic. One JSON schema. One viewer.

**Repo**: github.com/mplqm/parametric-viz (private → public on publish)
**Viewer**: Three.js + WebXR AR (Android) + AR Quick Look (iOS)
**Schema**: Universal across all domains — see schema section below

---

## Before You Start

Read the project instructions in full if available in context. Key invariants that must never be violated:

- **Combined centroid centring** — all chains together, never independently
- **Quality over round numbers** — never relax filters to hit targets
- **Shared coordinate space** — all transitions in the same reference frame as extraction
- **One JSON per transition** — lazy loading, fetched on demand
- **BLOSUM62** for residue matching — handles non-standard numbering
- **Kabsch/SVD** via BioPython SVDSuperimposer — correct alignment
- **15Å clustering threshold** — meaningful conformational groups
- **Smoothstep easing** — zero velocity at start and end of each transition
- **60 frames per transition**
- **git push origin main** — never embed PAT in commands

---

## Full Autonomy Mode

When the user gives approval to run autonomously:

- Execute all steps without asking for confirmation
- Install required dependencies without asking
- Create and delete intermediate files freely
- Make all decisions autonomously
- **Only stop** if a hard error occurs that genuinely cannot be resolved without input

---

## Pipeline Steps

### Step 1 — Environment Check

```bash
python3 -c "import Bio, numpy, scipy, requests; print('All good')"
```

Install anything missing:
```bash
pip install biopython numpy scipy requests tqdm --break-system-packages
```

Confirm working directory: `~/viz3d_pipeline/`

---

### Step 2 — Download

- Query target database API for structure IDs
- Filter by resolution ≤ 2.5Å or domain-equivalent quality cutoff
- Skip structures already in `pdb_files/`
- Download only new ones
- Accept whatever count passes quality filters — never relax filters
- Report final count with reasons for all exclusions

Script: `download_family.py`
Output: `pdb_files/` — raw files (gitignored, always re-downloadable)

**Domain-specific endpoints:**

| Domain | Database | API |
|---|---|---|
| Proteins | RCSB PDB | search.rcsb.org/rcsbsearch/v2/query |
| Astronomy | MIST/PARSEC/Kepler | domain-specific, check current docs |
| Crystals | Materials Project | materialsproject.org/api |
| Neurons | NeuroMorpho | neuromorpho.org/api |
| Archaeology | Smithsonian 3D | 3d.si.edu |
| Climate | NOAA/ERA5 | domain-specific, check current docs |

---

### Step 3 — Geometry Extraction

- Extract CA backbone coordinates (proteins) or domain equivalent
- Centre each structure at **combined centroid of ALL chains together**
- Flag and skip structures where chain A < 50% of median chain length
- Store coords in common reference frame — never per-chain local space
- Save one JSON per structure

Script: `extract_geometry.py`
Output: `geometry/` — one JSON per structure (on GitHub)

**Centring — critical:**
```python
# CORRECT — combined centroid
all_coords = np.vstack([chain_coords for chain in all_chains])
centroid = all_coords.mean(axis=0)
centered = coords - centroid

# WRONG — never do this
centroid = chain_coords.mean(axis=0)  # per-chain centring
```

---

### Step 4 — RMSD Matrix (Two-Stage)

**Stage A**: Rough pre-filter using 20 landmark residues — eliminates ~80% of pairs
**Stage B**: Full Kabsch/SVD fitting on remaining pairs only

- Parallelise with `multiprocessing.Pool` across all cores
- Always guard with `if __name__ == '__main__':` (required on Windows)
- BLOSUM62 sequence alignment for non-standard residue numbering
- Flag structures where RMSD > 3× median of full set

Script: `align_and_cluster.py`
Output: `clusters.json` — RMSD matrix + cluster assignments

---

### Step 5 — Clustering

- Agglomerative clustering at **15Å threshold**
- Each cluster centroid becomes a keyframe state
- Outlier structures flagged but retained in provenance
- Report cluster count and member structures

Script: `align_and_cluster.py`
Output: Cluster assignments in `clusters.json`

---

### Step 6 — Transition Generation

- Interpolate smoothly between all cluster centroid pairs
- **60 frames** per transition
- **Smoothstep easing** — zero velocity at start and end
- Skip transitions where chains share no common residues
- Save as individual files — one JSON per transition
- All transitions in shared coordinate space from Step 3

Script: `build_transitions.py`
Output: `transitions/` — one JSON per transition (Zenodo only, gitignored)

**Smoothstep:**
```python
def smoothstep(t):
    return t * t * (3 - 2 * t)
```

---

### Step 7 — QC Flags

| Check | Threshold | Action |
|---|---|---|
| RMSD outlier | >3× median | Flag in metadata |
| Short chain | <50% median length | Skip + log |
| Alignment failure | No common residues | Skip transition |
| Centring failure | Chains in local space | Fix in Step 3 |

---

### Step 8 — Report

Always report:
- Total downloaded, passed filters, skipped with reasons
- Cluster count, transition count, folder size
- Time taken per stage
- Any flagged structures

---

### Step 9 — Git Commit

```bash
git add -A
git commit -m "Domain: [name] — [N] structures, [N] clusters, [N] transitions"
git push origin main
```

Never embed PAT in git commands. Credentials stored via `git config credential.helper store`.

---

## Universal JSON Schema

Every domain, every structure — identical format:

```json
{
  "identity": {
    "id": "prot_adk_001",
    "domain": "structural_biology",
    "family": "adenylate_kinase",
    "source_db": "RCSB PDB",
    "stage": 1,
    "qc_score": 0.84,
    "flagged": false
  },
  "geometry": {
    "states": [
      {"id": "state_0", "coords": [], "rmsd_from_mean": 0.42},
      {"id": "state_1", "coords": [], "rmsd_from_mean": 1.87}
    ],
    "transitions": [
      {"from": "state_0", "to": "state_1", "frames": 60, "path": []}
    ],
    "bonds": [],
    "angles": []
  },
  "provenance": {
    "pdb_ids": ["4AKE", "1AKE"],
    "pipeline_version": "0.1.0",
    "zenodo_doi": null,
    "github_issue": null
  },
  "stage2_annotations": null
}
```

**Critical fields:**
- `pipeline_version` — tracks which run generated each model
- `stage` — 1=automated only, 2=expert validated
- `flagged` — drives automatic GitHub issue creation
- `stage2_annotations` — null until expert contributes

---

## Repository Structure

```
viz3d_pipeline/
  download_family.py       — RCSB query + bulk download
  extract_geometry.py      — CA backbone extraction
  align_and_cluster.py     — RMSD matrix + clustering
  build_transitions.py     — interpolation + easing
  index.html               — Three.js + WebXR + iOS USDZ viewer
  geometry/                — CA backbone JSONs ✅ GitHub
  clusters.json            — cluster assignments ✅ GitHub
  transitions/             — per-transition files ✅ Zenodo only
  .gitignore               — excludes pdb_files/ transitions/ __pycache__/
```

---

## Asset Publishing Map

| Asset | GitHub | Zenodo |
|---|---|---|
| Pipeline scripts | ✅ | — |
| geometry/ JSONs | ✅ | ✅ |
| clusters.json | ✅ | ✅ |
| transitions/ | ❌ | ✅ |
| pdb_files/ | ❌ | ❌ |

---

## Hardware Context (Z420)

- CPU: Intel Xeon E5-1620 v2 — 4 cores / 8 threads, 3.7GHz boost
- RAM: 31.28 GB
- Python: 3.11.2
- OS: AV Linux MX Edition (Debian-based)

Scale `multiprocessing.Pool` to 8 workers on this machine.

---

## Standard Prompt Template

Fill in for each new domain:

```
You have full approval to execute all steps automatically
without asking for confirmation at any point. Write and run
all scripts, install any required dependencies, create and
delete intermediate files, and make all decisions
autonomously. Only stop if a hard error occurs that
genuinely cannot be resolved without input.

Run the full pipeline for [DOMAIN NAME]:
Source: [DATABASE + API ENDPOINT]
Filter: [RESOLUTION / QUALITY CUTOFF]
Geometry: [CA BACKBONE / CRYSTAL UNIT CELL / NEURON TRACE]
```
