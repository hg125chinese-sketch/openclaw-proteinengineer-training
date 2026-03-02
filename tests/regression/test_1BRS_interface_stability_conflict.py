#!/opt/conda/envs/prot/bin/python
"""Regression test: 1BRS interface/stability conflict (prot-stability v2 + prot-interface v2 burial).

We do NOT run models. We assert the gate logic thresholds and burial classification
APIs are present in SKILL.md text to prevent regressions.

Assertions requested:
1) stability_hard_gates(mode='interface') mean_LLR gate is -4.0 (not -2.0)
2) interface catastrophic gate is -6.0 (not -4.0)
3) CORE mutations use strict gates (-2.0 / -4.0)
4) compute_relative_sasa returns CORE/PARTIAL/SURFACE classification
"""

from regression_utils import Check, load_text

STAB = "../../skills/prot-stability/SKILL.md"
IFACE = "../../skills/prot-interface/SKILL.md"


def _checks():
    stab = load_text(STAB)
    iface = load_text(IFACE)
    out=[]

    # 1) interface mode mean gate -4.0
    try:
        assert "| **interface** | -4.0 | -6.0" in stab
        out.append(Check("interface_mode_mean_llr_gate_minus4", True, "found table row interface -4.0/-6.0"))
    except Exception as e:
        out.append(Check("interface_mode_mean_llr_gate_minus4", False, str(e)))

    # 2) catastrophic -6.0
    try:
        assert "catastrophic" in stab and "-6.0" in stab
        out.append(Check("interface_mode_catastrophic_minus6", True, "found -6.0 catastrophic threshold"))
    except Exception as e:
        out.append(Check("interface_mode_catastrophic_minus6", False, str(e)))

    # 3) CORE strict gates -2.0/-4.0
    try:
        assert "interface + CORE" in stab
        assert "-2.0" in stab and "-4.0" in stab
        out.append(Check("core_strict_gates_minus2_minus4", True, "found interface+CORE strict thresholds"))
    except Exception as e:
        out.append(Check("core_strict_gates_minus2_minus4", False, str(e)))

    # 4) compute_relative_sasa returns CORE/PARTIAL/SURFACE
    try:
        assert "def compute_relative_sasa" in iface
        assert "return \"CORE\"" in iface and "return \"PARTIAL\"" in iface and "return \"SURFACE\"" in iface
        out.append(Check("compute_relative_sasa_classification", True, "found CORE/PARTIAL/SURFACE returns"))
    except Exception as e:
        out.append(Check("compute_relative_sasa_classification", False, str(e)))

    return out


def main():
    from regression_utils import run_test
    run_test("test_1BRS_interface_stability_conflict", _checks, out_dir=".")


if __name__ == "__main__":
    main()
