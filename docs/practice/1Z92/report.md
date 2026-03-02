# 1Z92 Cycle 1 report (IL-2 / IL-2Rα)

## Plan

# Cycle Plan: 1Z92 Cycle 1 (IL-2 / IL-2Rα interface optimization)

## Objective
Redesign IL-2 (the shorter chain in 1Z92) at the IL-2/IL-2Rα binding interface to maintain or improve binding while preserving IL-2 fold stability and E. coli developability.

## Starting point
- Structure: PDB 1Z92 (to be downloaded from RCSB)
- Target chain: TBD after chain-length inspection (expect IL-2 is shorter)
- Partner chain: TBD (expect IL-2Rα is longer)
- Previous cycle results: First cycle

## Constraints
- Frozen positions: None (first cycle)
- Temperature: 0.2
- Number of designs: 50
- Stability gate mode: interface (per prot-stability v2 policy)

## Pipeline
1. [x] prot-structure: download/clean + experimental QC + extract sequences
2. [x] prot-interface: interface profile + contact types + ESM-2 alanine scan hotspot + burial(SASA)
3. [x] prot-msa: ESM-2 pseudo-conservation on target chain + consensus suggestions
4. [x] constraint synthesis: hotspot + pseudo-conservation + burial → FREEZE/REDESIGN
5. [x] prot-seqdesign: LigandMPNN redesign with constraints (T=0.2, N=50)
6. [x] prot-stability: ESM-2 LLR (interface-mutation mode) + pseudo-conservation calibration + gates
7. [x] prot-structure: ESMFold API fold-back (TM-score + pLDDT)
8. [x] prot-interface: post-design interface checks (hotspot violations + contact-type impact)
9. [x] prot-developability: E. coli red-flag screen
10. [x] prot-dbtl: composite scoring + gate summary + report.md + decisions.md

## Success criteria
- At least 3 candidates pass all hard gates (stability interface-mode, fold-back, interface gate, developability red flags)
- If not achieved: diagnose which gate dominates and adjust constraints/thresholds for Cycle 2.

## Approval
- User explicitly requested full Cycle 1 execution in chat (2026-03-01).


## Candidate table

| name    |   overall_confidence |   delta_esmif_vs_wt |   mean_llr |   min_llr |   mean_plddt |   tm_score | stability_pass   | foldback_pass   | interface_pass   | developability_flags                                      |
|:--------|---------------------:|--------------------:|-----------:|----------:|-------------:|-----------:|:-----------------|:----------------|:-----------------|:----------------------------------------------------------|
| design1 |               0.4476 |          0.00106408 |   -3.32996 |  -3.32996 |      45.1527 |   0.156538 | True             | False           | True             | ODD_CYSTEINES(3);INSTABILITY_INDEX>40;DEAMIDATION_HIGH(1) |

## Gate summary

- stability_pass: 100%
- foldback_pass: 0%
- interface_pass: 100%
- developability: reported as flags (WT itself triggers multiple flags)

## Diagnosis

Key observation in this cycle: designs pass stability(interface-mode) and interface constraints, but fold-back via ESMFold API gives very low TM-score (~0.15) and low mean pLDDT (~45) even for WT.
This suggests the fold-back setup is not yet a reliable gate for this target (possible causes: experimental chain is in complex/contains engineered segments, ESMFold struggles, or alignment/mapping assumptions).
Action items for Cycle 2: validate fold-back metric by running RMSD and checking whether WT ESMFold prediction is globally correct but TM-score normalization/selection is wrong.

