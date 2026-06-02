# Contributing

This project publishes automated structural visualizations — Stage 1. The science inside them is real, but the meaning hasn't been assigned yet.

That's where you come in.

---

## What Stage 2 means

The pipeline extracts geometry from public databases and groups structures by shape. It doesn't know what those shapes mean. It can identify that a protein has two distinct conformations. It can't tell you that one is the ATP-bound state and the other is the apo form.

Domain experts can. That knowledge — cross-referenced against literature, validated against experimental data — is Stage 2.

Stage 2 contributions don't touch the geometry. They add meaning to it.

---

## What a contribution looks like

A single JSON block added to the `stage2_annotations` field of one or more structure files:

```json
"stage2_annotations": {
  "contributor": "your_github_username",
  "date": "2026-09-14",
  "cluster_labels": {
    "state_0": "open conformation — apo form",
    "state_1": "closed conformation — ATP-bound"
  },
  "notes": "State 1 corresponds to the closed intermediate described in Müller et al. 2019. State 2 shows partial closure consistent with AMP binding.",
  "references": ["doi:10.1038/s41586-019-XXXX-X"],
  "validated": true
}
```

That's it. No code required. No pipeline knowledge required. Just your domain expertise, written down in a structured way.

---

## How to contribute

### Option 1 — GitHub Pull Request

For researchers comfortable with Git:

1. Fork the repository at `github.com/mplqm/parametric-viz`
2. Find the structure file you want to annotate in `geometry/`
3. Add your `stage2_annotations` block
4. Open a pull request with a brief description of what you've annotated and why

The PR description doesn't need to be formal. A few sentences explaining what these conformational states represent is enough.

### Option 2 — GitHub Issue

If you'd rather not touch the files directly:

1. Open an issue at `github.com/mplqm/parametric-viz/issues`
2. Use the label `stage2-annotation`
3. Paste your annotation as a code block
4. Include which structure ID(s) you're annotating

We'll incorporate it and credit you in the file.

### Option 3 — Email

If neither of the above works:

Write your annotation in the format shown above and send it to the address in the repo profile. Same credit, same inclusion — just a slower process.

---

## What to annotate

### Cluster labels

Each structure file has a set of states — cluster centroids representing distinct conformational postures. Your job is to tell us what they are.

Look at the viewer. The states are the positions the structure holds between transitions. Open the structure on your phone or desktop, watch it move, and identify what you're seeing against your knowledge of this molecule.

Label each state with:
- A short name (`open conformation`, `ATP-bound closed form`)
- Optionally, a longer description in `notes`

### Notes

Free text. Anything useful — literature connections, caveats, known limitations of the automated clustering, structures you think were mis-assigned, interesting outliers.

### References

DOIs preferred. Any stable identifier works — PubMed ID, arXiv, book chapter. These get linked in the viewer so anyone exploring the structure can follow up.

---

## Flagged structures

Some structures are automatically flagged by the pipeline — RMSD outliers, unusual chain lengths, alignment failures. These have open GitHub issues.

Flagged structures are exactly the ones most likely to be scientifically interesting. If you recognise a flagged structure as a known unusual conformation, a disease-associated variant, or something worth keeping — say so in the issue. That context is valuable even if the automated QC marked it as anomalous.

---

## Credit

Every contribution is credited in the structure file itself — your GitHub username and the date. When the Zenodo dataset is updated to include Stage 2 annotations, contributors are listed in the dataset metadata.

The goal is to make it easy for domain communities to take ownership of their structures. The pipeline does the infrastructure work. The science belongs to the people who understand it.

---

## What we're not asking for

- You don't need to validate the geometry. The pipeline handles that.
- You don't need to run any code.
- You don't need to understand Three.js, WebXR, or how the viewer works.
- You don't need to annotate everything — partial contributions are welcome.
- You don't need institutional approval to contribute. This is open science.

---

## Questions

Open an issue with the label `question`. No question is too basic. The project is designed for non-expert users at the front end and expert contributors at the back — both matter equally.

---

*Built by one person and Claude, on a workstation from 2013.*
