# Parametric Viz — CLAUDE.md

## What This Project Is

A fully automated pipeline that transforms publicly available scientific structure data into interactive parametric 3D visualizations. Domain-agnostic — one JSON schema, one Three.js + WebXR viewer.

**Current state:** Proteins (adenylate kinase) are live. 942 structures → 313 clusters → 30,735 transitions. Viewer works on iOS (AR Quick Look) and Android (WebXR).

**Vision:** Same pipeline, same phone URL — proteins first, then stars, crystals, neurons, archaeological objects, climate models.

**Two-stage architecture:**
- Stage 1 — fully automated pipeline (download → geometry → RMSD → cluster → transitions → publish)
- Stage 2 — expert annotation layer (null until a domain expert contributes)

See `SKILL.md` for full pipeline reference and invariants. See `docs/` for schema, pipeline, publishing, and contributing docs.

---

## Session logging (required)

Call `log_session()` at the end of every session. Do **not** hand-write entries directly into the log file.

**Target file:** `/home/mp/viz3d_pipeline/docs/decisions/parametricviz_session_log.md`

**Format:** narrative style — `---` separator, `**timestamp — title**` heading, prose description, `**Next:** next_step` line.

**Import and call:**

```python
import sys
sys.path.insert(0, '/home/mp/viz3d_pipeline')
from session_log import log_session

log_session(
    title        = 'Short session title',
    description  = 'What was done and why — narrative prose.',
    next_step    = 'What to do next session.',
    current_file = 'the_main_file_worked_on.py',
    open_tasks   = [
        'Task one',
        'Task two',
    ],
    next_session = 'One sentence: what next session should start with.',
    notes        = '',   # optional
)
```
