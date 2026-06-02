# JSON Schema

Every structure, every domain, every pipeline run — identical format.

This schema is the contract between the pipeline and the viewer, between Stage 1 automation and Stage 2 expert annotation, and between this project and anyone who wants to build on top of it.

---

## Full Schema

```json
{
  "identity": {
    "id": "prot_adk_001",
    "domain": "structural_biology",
    "family": "adenylate_kinase",
    "source_db": "RCSB PDB",
    "stage": 1,
    "qc_score": 0.84,
    "flagged": false,
    "coord_units": "angstroms"
  },
  "geometry": {
    "states": [
      {
        "id": "state_0",
        "coords": [[x, y, z], ...],
        "rmsd_from_mean": 0.42
      },
      {
        "id": "state_1",
        "coords": [[x, y, z], ...],
        "rmsd_from_mean": 1.87
      }
    ],
    "transitions": [
      {
        "from": "state_0",
        "to": "state_1",
        "frames": 60,
        "path": "transitions/state_0_to_state_1.json"
      }
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

---

## Field Reference

### `identity`

| Field | Type | Description |
|---|---|---|
| `id` | string | Unique identifier. Format: `[domain_prefix]_[family]_[index]` |
| `domain` | string | Scientific domain. See domain values below. |
| `family` | string | Structural family within the domain. |
| `source_db` | string | Database this structure was downloaded from. |
| `stage` | integer | 1 = automated pipeline only. 2 = expert validated. |
| `qc_score` | float | Quality score 0.0–1.0. Derived from source database metrics. |
| `flagged` | boolean | True if any QC check failed. Triggers GitHub issue on publish. |

**Domain values:**

| Domain | `domain` string |
|---|---|
| Proteins | `structural_biology` |
| Astronomy | `astronomy` |
| Crystals | `materials_science` |
| Neurons | `neuroscience` |
| Archaeology | `archaeology` |
| Climate | `climate_science` |

**ID format examples:**
```
prot_adk_001       — protein, adenylate kinase, index 1
astr_solar_042     — astronomy, solar, index 42
crys_perovskite_7  — crystal, perovskite family, index 7
```

---

### `geometry`

#### `states`

Each state is a cluster centroid — a representative structure computed from a group of similar measured structures.

| Field | Type | Description |
|---|---|---|
| `id` | string | State identifier. Format: `state_[index]` |
| `coords` | array of [x, y, z] | 3D coordinates. One per residue (proteins) or domain equivalent. Centred at combined centroid of all chains. Units are domain-specific — see coordinate units below. |
| `rmsd_from_mean` | float | RMSD of this state from the mean structure of its cluster. Lower = more typical. |

**Coordinate frame:**
All coordinates are in the same reference frame — centred at the combined centroid of all chains together. This is established in `extract_geometry.py` and must never change between stages.

#### `transitions`

Each transition is a pointer to a separate file containing the interpolated frames between two states.

| Field | Type | Description |
|---|---|---|
| `from` | string | Source state ID |
| `to` | string | Destination state ID |
| `frames` | integer | Always 60. |
| `path` | string | Relative path to the transition JSON file. |

Transition files live in `transitions/` — Zenodo only, not on GitHub. The viewer fetches them on demand.

**Transition file format** (`transitions/state_0_to_state_1.json`):
```json
{
  "from": "state_0",
  "to": "state_1",
  "frames": [
    [[x, y, z], ...],
    [[x, y, z], ...],
    ...
  ]
}
```

60 entries in `frames`. Each entry is a full coordinate set — one [x, y, z] per residue. Smoothstep easing applied to the interpolation parameter.

#### `bonds`

Pairs of residue indices that share a bond. Used by the viewer to draw backbone tubes.

```json
"bonds": [[0, 1], [1, 2], [2, 3], ...]
```

For proteins, this is simply sequential CA pairs along the backbone. For other domains, bonds represent whatever structural connections are meaningful.

#### `angles`

Reserved. Empty array in Stage 1. May be populated by Stage 2 expert annotation.

---

### `provenance`

| Field | Type | Description |
|---|---|---|
| `pdb_ids` | array of strings | Source structure IDs from the originating database. |
| `pipeline_version` | string | Semver string of the pipeline run that generated this file. |
| `zenodo_doi` | string or null | DOI assigned after Zenodo upload. Null until published. |
| `github_issue` | string or null | GitHub issue URL if this structure was flagged. Null otherwise. |

`pipeline_version` is critical for reproducibility. If results change between runs, this field tells you which version produced which output.

---

### `stage2_annotations`

Null in all Stage 1 outputs. Populated by domain experts who contribute via GitHub pull request.

**When populated:**
```json
"stage2_annotations": {
  "contributor": "github_username",
  "date": "2026-09-14",
  "cluster_labels": {
    "state_0": "open conformation",
    "state_1": "closed conformation — ATP-bound"
  },
  "notes": "State 2 corresponds to the intermediate described in Smith et al. 2019.",
  "references": ["doi:10.1038/s41586-019-XXXX-X"],
  "validated": true
}
```

This field is the entire Stage 2 contribution surface. Experts annotate meaning without touching geometry. The pipeline geometry is never modified after Stage 1.

---

## Versioning

`pipeline_version` follows semver: `MAJOR.MINOR.PATCH`

| Change | Version bump |
|---|---|
| Schema field added or removed | MAJOR |
| New domain added | MINOR |
| Bug fix, parameter change | PATCH |

When `pipeline_version` changes, document what changed and why in the repo changelog. All files generated by a given version are internally consistent.

---


---

## Coordinate Units

Units are domain-specific and must be documented when a new domain is added:

| Domain | Units | Notes |
|---|---|---|
| Proteins | Ångströms (Å) | CA backbone, centred at combined centroid |
| Astronomy | TBD — normalised parameter space or physical units | Defined when domain is added |
| Crystals | Ångströms (Å) | Unit cell atomic positions |
| Neurons | Micrometres (µm) | Dendritic trace coordinates |
| Archaeology | Metres (m) or normalised | Depends on source scan resolution |
| Climate | Degrees lat/lon + altitude (m) | Or normalised parameter space |

Add a  field to  when adding a new domain to make units machine-readable:



## What changes between domains

**Nothing in the schema.** Every domain produces identical JSON structure. The coordinate semantics differ — proteins use Ångströms, astronomy may use parsecs or normalised parameter space — but the shape of the data is always the same.

Document the coordinate units and semantics for each domain in `docs/domains/[domain].md` as domains are added.

---

## What never changes

- The schema structure
- The coordinate frame (combined centroid centring)
- The transition file format (60 frames, smoothstep easing)
- The `stage2_annotations` field shape
- The `pipeline_version` field

These are the load-bearing walls. Everything else is configurable.

---

*Built by one person and Claude, on a workstation from 2013.*
