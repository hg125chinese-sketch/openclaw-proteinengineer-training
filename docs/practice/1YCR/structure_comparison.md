# 1YCR prot-structure: v2-like vs v3 comparison

Target: 1YCR_clean.pdb chain B (13 aa peptide)

## Key expected differences
- v3 scales API pLDDT to 0-100 (v2-like stays 0-1)
- v3 clash detection excludes bonded neighbors (min_seq_sep=3)
- v3 gate uses RMSD for short peptides (<30aa), not TM-score

## Table
| name     |   plddt_mean_v2 |   plddt_max_v2 | plddt_gate_v2   |   clashes_v2 |   tm_score_v2 |   rmsd_v2 | gate_v2   | primary_metric_v2   |   plddt_mean_v3 |   plddt_max_v3 | plddt_gate_v3   |   clashes_v3 |   tm_score_v3 |   rmsd_v3 | gate_v3   | primary_metric_v3   |
|:---------|----------------:|---------------:|:----------------|-------------:|--------------:|----------:|:----------|:--------------------|----------------:|---------------:|:----------------|-------------:|--------------:|----------:|:----------|:--------------------|
| wildtype |        0.776923 |           0.82 | REJECT          |           11 |      0.356585 |  1.21081  | FAIL      | TM-score            |         77.6923 |             82 | MODERATE        |            0 |      0.356585 |  1.21081  | PASS      | RMSD                |
| design1  |        0.840769 |           0.88 | REJECT          |           11 |      0.436355 |  0.79508  | FAIL      | TM-score            |         84.0769 |             88 | MODERATE        |            0 |      0.436355 |  0.79508  | PASS      | RMSD                |
| design2  |        0.830769 |           0.86 | REJECT          |           11 |      0.420721 |  0.770609 | FAIL      | TM-score            |         83.0769 |             86 | MODERATE        |            0 |      0.420721 |  0.770609 | PASS      | RMSD                |
| design3  |        0.849231 |           0.88 | REJECT          |           11 |      0.418668 |  0.72975  | FAIL      | TM-score            |         84.9231 |             88 | MODERATE        |            0 |      0.418668 |  0.72975  | PASS      | RMSD                |
