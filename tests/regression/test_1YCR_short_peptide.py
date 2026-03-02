#!/opt/conda/envs/prot/bin/python
"""Regression test: 1YCR short peptide edge cases (prot-structure v3/v4).

We do NOT call ESMFold API. We assert that prot-structure SKILL.md encodes
(and therefore future edits must preserve) the known short-peptide fixes.

Assertions requested:
1) ESMFold API pLDDT scaling to 0-100 (not 0-1)
2) Clash detection uses min_seq_sep=3 and short-peptide severe_threshold=3 (<3 clashes expectation)
3) Short peptides (<30aa) use RMSD as primary metric (not TM-score)
4) Short peptide pLDDT gate threshold uses >=70 (not >=75)
"""

from regression_utils import Check, load_text

SKILL = "../../skills/prot-structure/SKILL.md"


def _checks():
    t = load_text(SKILL)
    out = []

    # 1) pLDDT 0-1 -> 0-100 scaling
    try:
        assert "plddt.max() <= 1.0" in t and "plddt = plddt * 100.0" in t
        out.append(Check("plddt_scale_0to1_to_0to100", True, "found scaling logic"))
    except Exception as e:
        out.append(Check("plddt_scale_0to1_to_0to100", False, str(e)))

    # 2) detect_clashes min_seq_sep=3 and short-peptide severe_threshold=3
    try:
        assert "def detect_clashes" in t and "min_seq_sep=3" in t
        assert "severe_threshold = 3 if is_short_peptide else 10" in t
        out.append(Check("clash_min_seq_sep_3_and_peptide_threshold", True, "found min_seq_sep=3 and severe_threshold=3"))
    except Exception as e:
        out.append(Check("clash_min_seq_sep_3_and_peptide_threshold", False, str(e)))

    # 3) Short peptides use RMSD as primary metric
    try:
        assert "primary_metric" in t and '"RMSD" if is_short else "TM-score"' in t
        assert "Short peptides (<30 aa)" in t and "Uses CA RMSD as primary metric instead" in t
        out.append(Check("short_peptide_primary_metric_rmsd", True, "found RMSD primary metric language"))
    except Exception as e:
        out.append(Check("short_peptide_primary_metric_rmsd", False, str(e)))

    # 4) Short peptide pLDDT gate threshold >= 70
    try:
        assert "Short peptide fold-back OK" in t and "plddt >= 70" in t
        # also ensure standard protein PASS uses 75 so we are not confusing
        assert "tm >= 0.85 and plddt >= 75" in t
        out.append(Check("short_peptide_plddt_gate_ge_70", True, "found plddt>=70 for peptide and >=75 for proteins"))
    except Exception as e:
        out.append(Check("short_peptide_plddt_gate_ge_70", False, str(e)))

    return out


def main():
    # manual runner compatibility
    from regression_utils import run_test
    run_test("test_1YCR_short_peptide", _checks, out_dir=".")


if __name__ == "__main__":
    main()
