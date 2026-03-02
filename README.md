# PE Training v1 (ProteinEngineer)

This directory is a **sanitized, open-source-ready** snapshot of the ProteinEngineer (PE) training artifacts.

## What this is

PE training v1 produced:
- 7 protein-engineering skills (documented as SKILL.md files)
- 3 practice runs (1YCR, 1BRS, 1Z92) with written reports
- a regression test suite to prevent known edge-case fixes from regressing

## Contents

- `skills/` — the 7 skill documents:
  - prot-seqdesign
  - prot-structure
  - prot-stability
  - prot-interface
  - prot-msa
  - prot-developability
  - prot-dbtl

- `training-summary.md` — high-level training summary
- `docs/` — public-facing documentation:
  - `PLAYBOOK.en.md` — reproducibility/verifiability conventions (sanitized)
  - `skill-review.en.md` — what was validated + key outcomes
  - `practice/` — markdown reports from practice runs (no PDBs / large CSVs)

- `tests/regression/` — regression tests (no API calls; asserts gate logic invariants)

## What is intentionally NOT included

- Any vendor code (e.g., LigandMPNN) — use upstream:
  - https://github.com/dauparas/LigandMPNN
- Any PDB files or large CSVs from `runs/`
- Any machine-specific paths or personal identifiers

## How to run regression tests

From the repository root:

```bash
/opt/conda/envs/prot/bin/python tests/regression/run_all.py
```

A non-zero exit code indicates a regression.

## License

This directory contains documentation and lightweight test scripts.
If you publish it, add an explicit LICENSE at the repository root.
