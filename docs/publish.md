# Publishing

`publish.py` automates the full publish sequence — Zenodo upload, DOI retrieval, JSON patching, GitHub issues, and repo visibility. No manual steps required after running it.

---

## Prerequisites

### Tokens

Two environment variables required:

```bash
export ZENODO_TOKEN=your_zenodo_token
export GITHUB_TOKEN=your_github_pat
```

**Zenodo token** — create at zenodo.org → Account → Applications → Personal Access Tokens. Scopes required: `deposit:write`, `deposit:actions`.

**GitHub token** — create at github.com → Settings → Developer settings → Personal access tokens. Scopes required: `repo` (full), `public_repo`.

Store these in `~/.bashrc` or `~/.profile` so they persist across sessions. Never commit them to git.

### Domain config

Each domain needs a config file at `domain_config/<domain>.json`:

```json
{
  "title": "Adenylate Kinase — Conformational States and Transitions",
  "domain": "structural_biology",
  "family": "adenylate_kinase",
  "description": "Parametric 3D visualization of adenylate kinase conformational dynamics. 942 structures from RCSB PDB, quality-filtered at resolution ≤ 2.5Å. 313 conformational clusters, 30,735 smooth transitions. Stage 1 automated output — expert annotation welcome.",
  "keywords": ["protein structure", "adenylate kinase", "conformational dynamics", "AR visualization", "WebXR"],
  "creators": [
    {
      "name": "mplqm",
      "affiliation": "parametric-viz"
    },
    {
      "name": "Claude (Anthropic)",
      "affiliation": "Anthropic"
    }
  ],
  "pipeline_version": "0.1.0",
  "structure_count": 942,
  "cluster_count": 313,
  "transition_count": 30735
}
```

---

## Usage

### Always sandbox first

Zenodo sandbox is a full dry-run environment. Identical to production — same API, same sequence — but deposits are not public and no real DOI is issued.

```bash
python3 publish.py --domain proteins --sandbox
```

Check the output:
- Deposit created ✅
- Files uploaded ✅
- DOI returned (sandbox DOI, starts with `10.5072/`) ✅
- JSONs patched with DOI ✅
- Git commit and push ✅
- GitHub issues opened for flagged structures ✅

If any step fails, fix it before running production.

### Production run

```bash
python3 publish.py --domain proteins
```

This is irreversible. A published Zenodo deposit cannot be deleted — only a new version can supersede it.

### Flags

```bash
--domain    Required. Matches domain_config/<domain>.json
--sandbox   Use Zenodo sandbox instead of production
--skip-git  Skip git commit and push (useful if already committed)
--skip-issue  Skip GitHub issue creation
```

---

## What publish.py does — step by step

**Step 1 — Create Zenodo deposit**
Opens a new deposit draft via the Zenodo API. Returns a deposit ID and a file bucket URL.

**Step 2 — Upload files**
Uploads `transitions/` and `geometry/` to the deposit bucket. Large files are streamed in chunks. Progress shown per file.

**Step 3 — Set metadata**
Populates the deposit record from `domain_config/<domain>.json` — title, description, keywords, creators, pipeline version.

**Step 4 — Publish and retrieve DOI**
Triggers the publish action. Zenodo assigns a DOI immediately. The script captures it.

**Step 5 — Patch JSONs**
Updates `zenodo_doi` in every JSON file in `geometry/`. Uses the DOI retrieved in Step 4.

**Step 6 — Git commit and push**
Commits the DOI-patched JSONs with message: `Publish: [domain] — DOI [doi]`. Pushes to origin main.

**Step 7 — Open GitHub issues**
Finds all geometry JSONs where `flagged: true`. Opens one GitHub issue per flagged structure using the Stage 2 invitation template. Labels: `stage-2`, `annotation-request`, `[domain]`.

**Step 8 — Make repo public**
Updates the GitHub repo visibility to public via the GitHub API.

---

## After publishing

Update the project instructions with:
- Zenodo DOI
- Date published
- Structure count, cluster count, transition count

Share the viewer URL with at least one domain expert community before wider announcement. This gives the community a chance to engage with Stage 2 before the project gets broader attention.

---

## Versioning

Each domain publish is a separate Zenodo deposit with its own DOI. If you re-run the pipeline with updated data, create a new Zenodo version rather than a new deposit — this preserves the DOI history and keeps citations stable.

```bash
python3 publish.py --domain proteins --new-version
```

(Flag to be implemented when first update is needed.)

---

## Troubleshooting

**Upload stalls on large files** — transitions/ can be 6GB+. Run overnight. The script resumes from the last successful file if interrupted (re-run the same command).

**DOI patch fails** — check write permissions on geometry/ files. Re-run with `--skip-git` to retry just the patch and commit without re-uploading.

**GitHub issue creation fails** — check GITHUB_TOKEN scopes. `repo` scope required. Run with `--skip-issue` to complete publish without issues, then open them manually.

**Repo visibility change fails** — GitHub API requires the token owner to be a repo admin. Confirm token permissions and retry.

---

*Built by one person and Claude, on a workstation from 2013.*
