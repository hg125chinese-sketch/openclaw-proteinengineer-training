# 1YCR interface-only sequence design report

## Inputs
- PDB: 1YCR (cleaned chains A=MDM2, B=p53 peptide)
- Interface residues designed (B within 5Å of A): B17 B18 B19 B20 B22 B23 B25 B26 B27 B28 B29
- LigandMPNN: protein_mpnn, temperature=0.1, 50 samples

## Wildtype (from LigandMPNN FASTA first record)
- Chain B WT sequence: `ETFSDLWKLLPEN`
- WT overall_confidence: 0.4201 (NOTE: WT header did not include it; this is a fallback median if missing)
- WT ESM-IF ll_fullseq (chain B): -2.3069
- WT ESM-IF ll_withcoord (chain B): -2.3069

## LigandMPNN designs summary (n=50)
- overall_confidence: mean=0.4215, best=0.4404
- interface recovery frac (computed): mean=0.631, min=0.545, max=0.636
- interface mutation count: mean=4.06, min=4, max=5

## Hard gates (Phase 4.1)
- Passed: 25/50 (50.0%)

## ESM-IF cross-validation on Top 5 (Phase 3)
|   sample_idx |   overall_confidence |   iface_mut_count |   iface_recovery_frac |   esmif_ll_fullseq |   delta_esmif_ll_fullseq_vs_wt |   esmif_ll_withcoord |   delta_esmif_ll_withcoord_vs_wt | agreement   |
|-------------:|---------------------:|------------------:|----------------------:|-------------------:|-------------------------------:|---------------------:|---------------------------------:|:------------|
|            9 |               0.4404 |                 4 |              0.636364 |           -1.96848 |                       0.338422 |             -1.96848 |                         0.338422 | BOTH_BETTER |
|           50 |               0.4372 |                 4 |              0.636364 |           -1.96848 |                       0.338422 |             -1.96848 |                         0.338422 | BOTH_BETTER |
|           12 |               0.4371 |                 4 |              0.636364 |           -1.96848 |                       0.338422 |             -1.96848 |                         0.338422 | BOTH_BETTER |
|           30 |               0.4341 |                 4 |              0.636364 |           -1.96848 |                       0.338422 |             -1.96848 |                         0.338422 | BOTH_BETTER |
|           36 |               0.432  |                 4 |              0.636364 |           -1.96848 |                       0.338422 |             -1.96848 |                         0.338422 | BOTH_BETTER |

## Notes
- This LigandMPNN build outputs `overall_confidence` + `seq_rec` instead of `score/global_score` described in SKILL.md.
- Hard gates and ranking were adapted accordingly; see skill-feedback.md for suggested updates.
