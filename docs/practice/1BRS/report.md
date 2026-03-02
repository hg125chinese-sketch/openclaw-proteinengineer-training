# 1BRS integration: prot-structure → prot-interface → prot-seqdesign → prot-stability

## Step 1 (prot-structure): experimental structure QC
- PDB: 1BRS (chains A=barnase, D=barstar)
- Chain lengths: A=108, D=87
- Missing backbone atoms: A=0, D=0
- Clash check (min_seq_sep=3): n_clashes=0

## Step 2 (prot-interface): interface + hotspots + constraints
- Chain D interface residues: 19
- STRONG_HOTSPOT: 17 (89.5%)  (expected ~10–30% typically)
- Constraints: FREEZE=17, REDESIGN=2
  - REDESIGN residues: D36 D40

## Step 3 (prot-seqdesign): LigandMPNN redesign (T=0.2, N=50) + ESM-IF
- Designed residues: D36, D40
- Unique sequences in Top5: 1
- WT ESM-IF avg_ll: -1.3199

Top5 (by MPNN confidence / ESM-IF delta):
|   overall_confidence |   design_recovery |   delta_esmif_vs_wt |
|---------------------:|------------------:|--------------------:|
|               0.6186 |               0.5 |          -0.0182275 |
|               0.6121 |               0.5 |          -0.0182276 |
|               0.6204 |               0.5 |          -0.0182277 |
|               0.614  |               0.5 |          -0.0182277 |
|               0.6155 |               0.5 |          -0.0182279 |

## Step 4 (prot-stability): ESM-2 LLR + consensus + hard gates
| name   |   n_mutations |   mean_llr |   min_llr |   consensus_compat_rate | gate_pass   | gate_fail_reasons                           |
|:-------|--------------:|-----------:|----------:|------------------------:|:------------|:--------------------------------------------|
| cand1  |             1 |   -6.19224 |  -6.19224 |                       1 | False       | mean_LLR<-2.0;min_LLR<-4.0;>70% unfavorable |

## Step 5 (prot-structure, forced): ESMFold API fold-back (Top5 even if stability FAIL)
| name   |   pred_mean_plddt |   tm_score_vs_exp_chainD |   pred_len_ca |   exp_len_ca |
|:-------|------------------:|-------------------------:|--------------:|-------------:|
| cand1  |           87.2068 |                 0.970365 |            87 |           87 |

## Step 6 (prot-interface, forced): hotspot constraint gate (Top5 even if stability FAIL)
| name   |   n_mutations | passed   | mutations   | violations   | warnings   |
|:-------|--------------:|:---------|:------------|:-------------|:-----------|
| cand1  |             1 | True     | A36S        |              |            |

## Step 7: Summary table (diagnostic)
| name   |   mpnn_overall_confidence |   esmif_avg_ll |   delta_esmif_vs_wt |   mean_llr |   min_llr |   pred_mean_plddt |   tm_score_vs_exp_chainD | interface_gate_pass   | interface_mutations   |
|:-------|--------------------------:|---------------:|--------------------:|-----------:|----------:|------------------:|-------------------------:|:----------------------|:----------------------|
| cand1  |                    0.6186 |       -1.33816 |          -0.0182275 |   -6.19224 |  -6.19224 |           87.2068 |                 0.970365 | True                  | A36S                  |

## Diagnosis: why stability FAIL?
In this run, the candidate(s) show **fold-back PASS** (high pLDDT, TM-score~1) and **interface gate PASS**, but **ESM-2 LLR is catastrophically negative**.
This pattern is consistent with ESM-2 LLR acting as an **evolutionary plausibility prior** (sequence likelihood under a broad MSA-trained model), not a pure physical stability metric.
Actionable next tests:
- Verify whether the mutated position is buried/core vs solvent/interface (LLR is often harsher for core).
- Consider relaxing LLR hard gates for interface-only mutations, or switching to a stability proxy calibrated on barstar family (MSA-based or physics-based).

## Final recommendation
No candidates pass the current prot-stability hard gates (ESM-2 LLR).
However, fold-back + interface checks suggest the design is structurally plausible; the blocker is the LLR gate.
