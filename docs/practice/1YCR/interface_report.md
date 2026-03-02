# prot-interface practice: 1YCR (p53-MDM2)

PDB: <ABS_PATH_REDACTED> (A=MDM2, B=p53 peptide) cutoff=5.0Å

## Phase 1: Interface profile (chain B)
|   resseq_b | resname_b   |   n_contacts |   n_hbond |   n_salt |   n_aromatic |   n_hydrophobic |   min_distance |
|-----------:|:------------|-------------:|----------:|---------:|-------------:|----------------:|---------------:|
|         19 | PHE         |           88 |         1 |        0 |           16 |              28 |        3.01698 |
|         29 | ASN         |           84 |         3 |        0 |            0 |              12 |        2.6809  |
|         23 | TRP         |           68 |         1 |        0 |           17 |              22 |        2.83101 |
|         28 | GLU         |           48 |         1 |        0 |            1 |              10 |        3.16303 |
|         17 | GLU         |           41 |         0 |        1 |            0 |               3 |        3.10687 |
|         26 | LEU         |           32 |         0 |        0 |            9 |               7 |        3.36003 |
|         27 | PRO         |           31 |         0 |        0 |            3 |              12 |        3.24413 |
|         22 | LEU         |           28 |         0 |        0 |            3 |              14 |        3.3293  |
|         18 | THR         |           19 |         0 |        0 |            0 |               7 |        3.59936 |
|         25 | LEU         |           13 |         1 |        0 |            0 |               2 |        3.16947 |
|         20 | SER         |            9 |         1 |        0 |            0 |               2 |        3.09368 |

## Phase 2: Alanine-scan hotspots (chain B interface positions)
|   position | wt_aa   |     wt_ll |   ala_ll |   ala_llr | category         |
|-----------:|:--------|----------:|---------:|----------:|:-----------------|
|         10 | L       | -0.209393 | -4.64941 |  -4.44002 | STRONG_HOTSPOT   |
|          6 | L       | -0.219408 | -4.45355 |  -4.23415 | STRONG_HOTSPOT   |
|          9 | L       | -0.278291 | -4.35464 |  -4.07635 | STRONG_HOTSPOT   |
|          4 | S       | -0.337548 | -4.06907 |  -3.73152 | STRONG_HOTSPOT   |
|         12 | E       | -0.445982 | -4.01736 |  -3.57138 | STRONG_HOTSPOT   |
|          1 | E       | -1.24189  | -4.56385 |  -3.32197 | STRONG_HOTSPOT   |
|          2 | T       | -0.473983 | -3.78157 |  -3.30758 | STRONG_HOTSPOT   |
|         13 | N       | -0.53382  | -3.74635 |  -3.21253 | STRONG_HOTSPOT   |
|          3 | F       | -0.514258 | -3.68814 |  -3.17388 | STRONG_HOTSPOT   |
|         11 | P       | -0.697554 | -3.7408  |  -3.04325 | STRONG_HOTSPOT   |
|          7 | W       | -1.21452  | -3.29257 |  -2.07805 | MODERATE_HOTSPOT |

## Phase 2+1 merged: Hotspot + contact profile
|   position | wt_aa   |   ala_llr | category         |   n_contacts |   n_hbond |   n_salt |   n_hydrophobic |   min_distance |
|-----------:|:--------|----------:|:-----------------|-------------:|----------:|---------:|----------------:|---------------:|
|         10 | L       |  -4.44002 | STRONG_HOTSPOT   |           32 |         0 |        0 |               7 |        3.36003 |
|          6 | L       |  -4.23415 | STRONG_HOTSPOT   |           28 |         0 |        0 |              14 |        3.3293  |
|          9 | L       |  -4.07635 | STRONG_HOTSPOT   |           13 |         1 |        0 |               2 |        3.16947 |
|          4 | S       |  -3.73152 | STRONG_HOTSPOT   |            9 |         1 |        0 |               2 |        3.09368 |
|         12 | E       |  -3.57138 | STRONG_HOTSPOT   |           48 |         1 |        0 |              10 |        3.16303 |
|          1 | E       |  -3.32197 | STRONG_HOTSPOT   |           41 |         0 |        1 |               3 |        3.10687 |
|          2 | T       |  -3.30758 | STRONG_HOTSPOT   |           19 |         0 |        0 |               7 |        3.59936 |
|         13 | N       |  -3.21253 | STRONG_HOTSPOT   |           84 |         3 |        0 |              12 |        2.6809  |
|          3 | F       |  -3.17388 | STRONG_HOTSPOT   |           88 |         1 |        0 |              28 |        3.01698 |
|         11 | P       |  -3.04325 | STRONG_HOTSPOT   |           31 |         0 |        0 |              12 |        3.24413 |
|          7 | W       |  -2.07805 | MODERATE_HOTSPOT |           68 |         1 |        0 |              22 |        2.83101 |

## Phase 3: Design constraints
- FREEZE: B26 B22 B25 B20 B28 B17 B18 B29 B19 B27 B23
- REDESIGN: 
- CAUTION: 

## Phase 6: Interface gates on prior Top3 designs
```json
{
  "design1": {
    "mpnn_overall_confidence": 0.4404,
    "mutations": [
      "S20E",
      "L25K",
      "E28Q",
      "N29S"
    ],
    "passed": false,
    "violations": [
      "S20E mutates STRONG_HOTSPOT (ala_LLR=-3.73)",
      "L25K mutates STRONG_HOTSPOT (ala_LLR=-4.08)",
      "E28Q mutates STRONG_HOTSPOT (ala_LLR=-3.57)",
      "N29S mutates STRONG_HOTSPOT (ala_LLR=-3.21)"
    ],
    "warnings": []
  },
  "design2": {
    "mpnn_overall_confidence": 0.4156,
    "mutations": [
      "T18S",
      "S20E",
      "L25K",
      "E28Q",
      "N29S"
    ],
    "passed": false,
    "violations": [
      "T18S mutates STRONG_HOTSPOT (ala_LLR=-3.31)",
      "S20E mutates STRONG_HOTSPOT (ala_LLR=-3.73)",
      "L25K mutates STRONG_HOTSPOT (ala_LLR=-4.08)",
      "E28Q mutates STRONG_HOTSPOT (ala_LLR=-3.57)",
      "N29S mutates STRONG_HOTSPOT (ala_LLR=-3.21)"
    ],
    "warnings": []
  },
  "design3": {
    "mpnn_overall_confidence": 0.3987,
    "mutations": [
      "S20E",
      "L25Q",
      "E28Q",
      "N29S"
    ],
    "passed": false,
    "violations": [
      "S20E mutates STRONG_HOTSPOT (ala_LLR=-3.73)",
      "L25Q mutates STRONG_HOTSPOT (ala_LLR=-4.08)",
      "E28Q mutates STRONG_HOTSPOT (ala_LLR=-3.57)",
      "N29S mutates STRONG_HOTSPOT (ala_LLR=-3.21)"
    ],
    "warnings": []
  }
}
```
