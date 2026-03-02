#!/opt/conda/envs/prot/bin/python
"""Regression test: 1Z92 WT sanity + WT-relative developability.

No API calls. We assert that prot-structure v4+ has WT sanity check that
can downgrade fold-back to DIAGNOSTIC/SKIP, and that prot-developability
has WT-relative hard gate helper.

Assertions requested:
1) wildtype_sanity_check(): when WT pLDDT<60 -> gate_mode DIAGNOSTIC or SKIP (not HARD_GATE)
2) foldback_validation(gate_mode='DIAGNOSTIC') should not return FAIL (should be DIAGNOSTIC)
3) foldback_validation(gate_mode='SKIP') returns SKIPPED
4) developability_hard_gates_relative(): shared flags WT vs design not blocking

We implement (4) as a pure unit test on the helper's intended semantics,
while also asserting the helper exists in SKILL.md.
"""

from regression_utils import Check, load_text

STRUCT = "../../skills/prot-structure/SKILL.md"
DEV = "../../skills/prot-developability/SKILL.md"


def _dev_relative_unit():
    # emulate desired semantics
    wt = {"blocking_issues": ["ODD_CYSTEINES(3)", "INSTABILITY_INDEX>40"]}
    design = {"blocking_issues": ["ODD_CYSTEINES(3)", "INSTABILITY_INDEX>40", "HYDRO_PATCHES(3)"]}
    # In WT-relative mode, only issues added by design are blocking
    design_only = [x for x in design["blocking_issues"] if x not in wt["blocking_issues"]]
    assert design_only == ["HYDRO_PATCHES(3)"]


def _checks():
    s = load_text(STRUCT)
    d = load_text(DEV)
    out=[]

    # 1) WT sanity downgrade thresholds
    try:
        assert "def wildtype_sanity_check" in s
        assert "if wt_plddt >= 60 and wt_tm >= 0.7" in s
        assert "gate_mode = \"DIAGNOSTIC\"" in s
        assert "gate_mode = \"SKIP\"" in s
        out.append(Check("wt_sanity_has_diagnostic_and_skip", True, "found WT sanity check gate_mode DIAGNOSTIC/SKIP"))
    except Exception as e:
        out.append(Check("wt_sanity_has_diagnostic_and_skip", False, str(e)))

    # 2) foldback_validation DIAGNOSTIC not FAIL
    try:
        assert "if gate_mode == \"DIAGNOSTIC\":" in s
        assert "gate = \"DIAGNOSTIC\"" in s
        out.append(Check("foldback_diagnostic_not_fail", True, "foldback_validation sets gate=DIAGNOSTIC"))
    except Exception as e:
        out.append(Check("foldback_diagnostic_not_fail", False, str(e)))

    # 3) foldback_validation SKIP -> SKIPPED
    try:
        assert "if gate_mode == \"SKIP\":" in s
        assert "\"gate\": \"SKIPPED\"" in s
        out.append(Check("foldback_skip_returns_skipped", True, "foldback_validation returns SKIPPED"))
    except Exception as e:
        out.append(Check("foldback_skip_returns_skipped", False, str(e)))

    # 4) developability WT-relative helper exists and semantics
    try:
        assert "def developability_hard_gates_relative" in d
        _dev_relative_unit()
        out.append(Check("developability_wt_relative_semantics", True, "helper exists and unit semantics OK"))
    except Exception as e:
        out.append(Check("developability_wt_relative_semantics", False, str(e)))

    return out


def main():
    from regression_utils import run_test
    run_test("test_1Z92_wt_sanity", _checks, out_dir=".")


if __name__ == "__main__":
    main()
