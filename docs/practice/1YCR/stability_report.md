# prot-stability practice: 1YCR chain B (Top3 designs)

WT: `ETFSDLWKLLPEN`

## Summary (per design)
| name    |   mpnn_overall_confidence |   n_mutations |   mean_llr |   min_llr |   consensus_compat_rate | gate_pass   | gate_fail_reasons                           |
|:--------|--------------------------:|--------------:|-----------:|----------:|------------------------:|:------------|:--------------------------------------------|
| design2 |                    0.4156 |             5 |   -3.28366 |  -3.60355 |                    0.6  | False       | mean_LLR<-2.0;>70% unfavorable              |
| design1 |                    0.4404 |             4 |   -3.34086 |  -3.60355 |                    0.5  | False       | mean_LLR<-2.0;>70% unfavorable              |
| design3 |                    0.3987 |             4 |   -3.51292 |  -4.21347 |                    0.25 | False       | mean_LLR<-2.0;min_LLR<-4.0;>70% unfavorable |

## Worst single mutation (by LLR)
- design1: E12Q (LLR=-3.604)
- design2: E12Q (LLR=-3.604)
- design3: L9Q (LLR=-4.213)

## Ranking comparison
- prot-seqdesign MPNN (overall_confidence) rank: ['design1', 'design2', 'design3']
- prot-stability ESM-2 mean_LLR rank: ['design2', 'design1', 'design3']
- Agreement?: False
