# ProteinEngineer Skill Review (Public, English)

This document summarizes what was trained/validated in **PE training v1** and what the integration runs demonstrated.

## Skills included (7)

- **prot-seqdesign**: LigandMPNN-driven sequence design with ESM-IF cross-validation.
- **prot-structure**: ESMFold (API-first) structure scoring + fold-back validation, including short-peptide handling and WT sanity check.
- **prot-stability**: ESM-2 LLR screening, including an interface-mutation mode and conflict handling.
- **prot-interface**: Interface extraction, contact typing, alanine-scan hotspot proxy, and constraint generation.
- **prot-msa**: ESM-2 pseudo-conservation and consensus-style guidance to calibrate risk.
- **prot-developability**: Sequence-based red-flag screening for experimental risk.
- **prot-dbtl**: Cycle directory standard + orchestration + candidate tables and decisions.

## Practice runs included

### 1YCR (p53–MDM2) — short peptide edge cases

Purpose: validate the pipeline on a 13-aa peptide binder.

Key outcomes:
- Short peptide fold-back must use **RMSD** as primary metric (TM-score is unstable for <30 aa).
- ESMFold API pLDDT stored in B-factors may arrive in a 0–1 scale and must be rescaled to 0–100.
- Clash detection must exclude bonded neighbors (`min_seq_sep=3`) to avoid inflated clash counts.

### 1BRS (barnase–barstar) — interface/stability conflict diagnostic

Purpose: validate end-to-end interface redesign logic and detect proxy conflicts.

Key outcomes:
- ESM-2 LLR can strongly penalize interface mutations that otherwise look structurally plausible.
- The interface-mutation mode in stability is needed to avoid over-rejecting candidates.
- Burial classification (CORE/PARTIAL/SURFACE) is required to decide when strict gates should apply.

### 1Z92 (IL-2 / IL-2Rα) — full Cycle 1 (7-skill integration)

Purpose: run a prot-dbtl Cycle 1 directory-standard pipeline end-to-end.

Key outcomes:
- Hotspot+burial+conservation constraints can become too restrictive (few REDESIGN sites) and collapse sequence diversity.
- **WT fold-back sanity** is mandatory: if ESMFold cannot fold WT well for a target (low pLDDT/TM), fold-back must be downgraded to **diagnostic** for that target/cycle.
- Developability signals should often be interpreted **relative to WT baseline** for proteins whose WT already carries red flags.

## Regression tests

This release includes a small regression suite under `tests/regression/` to prevent known fixes from regressing:

- 1YCR short-peptide fold-back logic (pLDDT scaling, RMSD primary metric, peptide thresholds)
- 1BRS interface-mutation stability thresholds and burial classification APIs
- 1Z92 WT sanity check and WT-relative developability helper presence/semantics
