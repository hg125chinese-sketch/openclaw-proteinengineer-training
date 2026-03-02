---
name: prot-developability
description: Developability assessment for protein designs — expression risk, solubility, aggregation propensity, chemical liabilities (deamidation, oxidation, isomerization), net charge/pI, hydrophobic patches, free cysteines, and sequence complexity. Produces a "red-light table" of issues that could cause experimental failure. Triggers whenever designs are about to go to experimental validation and you need to flag manufacturability risks.
homepage: https://github.com/facebookresearch/esm
metadata: { "openclaw": { "emoji": "🧪", "requires": { "bins": ["python3"], "python": ["numpy", "pandas", "biopython"] } } }
---

# Protein Developability: The Last Gate Before the Bench

Your design passed fold-back, stability scoring, and interface gates. But the most common reason proteins fail in the lab is not bad design — it's bad developability: low expression, aggregation, chemical degradation, or insolubility. This skill catches those problems before you waste time and reagents.

## When to Use

- **Before ordering DNA** — screen final candidates for red flags
- After prot-dbtl composite scoring, as the **final filter**
- When choosing between candidates with similar computational scores
- When a design **expressed poorly** in a previous cycle and you need to understand why
- When designing for a specific **expression system** (E. coli vs CHO vs cell-free)

## Core Philosophy

1. **Developability is a veto, not a score.** A single red flag (e.g., unpaired cysteine, 5 adjacent Asn-Gly sites) can kill a design regardless of how good its other metrics are. Treat this as a pass/fail checklist, not a ranking tool.
2. **Sequence-based screens are fast and catch 80% of problems.** You don't need a crystal structure to flag deamidation sites or count hydrophobic patches. Run these checks on every candidate.
3. **Expression system matters.** A free cysteine is fine in E. coli cytoplasm (reducing) but lethal in CHO secretion (oxidizing). Always specify the intended expression system.
4. **Fix what you can, flag what you can't.** Some liabilities can be engineered out (Asn→Gln to remove deamidation). Others require accepting the risk (a critical catalytic Cys). Document both.
5. **This is the cheapest gate.** Every check here runs in milliseconds on CPU. There's no excuse for skipping it.

## Phase 0: Environment Verification

```python
#!/opt/conda/envs/prot/bin/python
"""Phase 0: Check dependencies."""
import importlib

print("=== prot-developability Phase 0 ===\n")

deps = {"numpy": None, "pandas": None, "Bio.SeqUtils.ProtParam": "biopython"}
missing = []
for mod, pkg in deps.items():
    try:
        importlib.import_module(mod.split(".")[0])
    except ImportError:
        missing.append(pkg or mod)

if missing:
    print(f"❌ Missing: {missing}")
else:
    print("✅ All dependencies OK (CPU only, no GPU needed)")

print("\n=== Phase 0 complete ===")
```

## Phase 1: Sequence Properties

### 1.1 Basic Properties (Length, MW, pI, Charge, GRAVY)

```python
#!/opt/conda/envs/prot/bin/python
"""Compute basic sequence properties."""
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

def sequence_properties(sequence):
    """Compute basic biophysical properties from sequence.

    Args:
        sequence: amino acid string (single-letter, no gaps)

    Returns:
        dict with MW, pI, charge_at_pH7, gravy, instability_index, etc.
    """
    # Clean sequence
    seq = sequence.upper().replace("X", "").replace("-", "").replace(".", "")
    pa = ProteinAnalysis(seq)

    mw = pa.molecular_weight()
    pi = pa.isoelectric_point()
    gravy = pa.gravy()  # Grand Average of Hydropathy (positive = hydrophobic)
    instability = pa.instability_index()
    aromaticity = pa.aromaticity()

    # Charge at pH 7
    charge_7 = pa.charge_at_pH(7.0)

    # Amino acid composition
    aa_comp = pa.get_amino_acids_percent()

    # Fraction charged
    frac_charged = sum(aa_comp.get(aa, 0) for aa in "DEKRH")

    # Fraction hydrophobic
    frac_hydrophobic = sum(aa_comp.get(aa, 0) for aa in "AVILMFWP")

    result = {
        "length": len(seq),
        "mw_da": float(mw),
        "pI": float(pi),
        "charge_pH7": float(charge_7),
        "gravy": float(gravy),
        "instability_index": float(instability),
        "aromaticity": float(aromaticity),
        "frac_charged": float(frac_charged),
        "frac_hydrophobic": float(frac_hydrophobic),
    }

    # Flags
    flags = []
    if instability > 40:
        flags.append("UNSTABLE (instability_index > 40)")
    if abs(charge_7) < 1 and len(seq) > 50:
        flags.append(f"NEAR_NEUTRAL_CHARGE (charge={charge_7:.1f} at pH 7) — aggregation risk")
    if gravy > 0:
        flags.append(f"HYDROPHOBIC (GRAVY={gravy:.2f} > 0) — solubility risk")
    if pi > 6.5 and pi < 7.5:
        flags.append(f"pI_NEAR_NEUTRAL ({pi:.1f}) — precipitation risk at physiological pH")

    result["flags"] = flags
    return result
```

### 1.2 Amino Acid Composition Flags

```python
def composition_flags(sequence):
    """Flag unusual amino acid composition.

    Returns:
        list of flag strings
    """
    seq = sequence.upper()
    L = len(seq)
    flags = []

    # Count specific amino acids
    counts = {aa: seq.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}

    # Free cysteines
    n_cys = counts.get("C", 0)
    if n_cys % 2 != 0:
        flags.append(f"ODD_CYSTEINES ({n_cys}) — unpaired Cys risk")
    if n_cys > 0:
        flags.append(f"HAS_CYSTEINES ({n_cys}) — check disulfide pairing")

    # Methionine (oxidation-sensitive)
    n_met = counts.get("M", 0)
    if n_met > L * 0.05:
        flags.append(f"HIGH_MET ({n_met}/{L}, {100*n_met/L:.1f}%) — oxidation risk")

    # Tryptophan (oxidation-sensitive)
    n_trp = counts.get("W", 0)
    if n_trp > L * 0.04:
        flags.append(f"HIGH_TRP ({n_trp}/{L}) — oxidation/fluorescence interference")

    # Proline content (affects folding kinetics)
    n_pro = counts.get("P", 0)
    if n_pro > L * 0.10:
        flags.append(f"HIGH_PRO ({n_pro}/{L}, {100*n_pro/L:.1f}%) — may slow folding")

    # Low complexity (runs of same amino acid)
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        for run_len in range(5, 20):
            if aa * run_len in seq:
                flags.append(f"HOMOPOLYMER_RUN ({aa}×{run_len}) — expression/folding risk")
                break

    return flags
```

## Phase 2: Chemical Liabilities

### 2.1 Deamidation Sites (Asn-X)

```python
#!/opt/conda/envs/prot/bin/python
"""Identify chemical liability motifs in protein sequences."""
import re

# Deamidation risk: Asn followed by small/flexible residues
_DEAMIDATION_HIGH = {"NG", "NS", "NT", "ND", "NH"}  # highest risk
_DEAMIDATION_MODERATE = {"NA", "NQ", "NE", "NK", "NR"}  # moderate

def find_deamidation_sites(sequence):
    """Find Asn deamidation-prone sites.

    Asn deamidation (Asn→Asp/isoAsp) is the most common chemical
    degradation. Risk depends on the +1 residue:
    - NG, NS, NT, ND, NH: HIGH risk (small/flexible +1)
    - NA, NQ, NE, NK, NR: MODERATE risk
    - NP, NW, NF, NY, NI, NL, NV: LOW risk (bulky +1 blocks)

    Returns:
        list of dicts with position, motif, risk level
    """
    sites = []
    seq = sequence.upper()

    for i in range(len(seq) - 1):
        if seq[i] != "N":
            continue
        dipeptide = seq[i:i+2]
        if dipeptide in _DEAMIDATION_HIGH:
            sites.append({
                "position": i + 1,
                "motif": dipeptide,
                "risk": "HIGH",
                "note": f"Asn{i+1}-{seq[i+1]}{i+2}: fast deamidation",
            })
        elif dipeptide in _DEAMIDATION_MODERATE:
            sites.append({
                "position": i + 1,
                "motif": dipeptide,
                "risk": "MODERATE",
                "note": f"Asn{i+1}-{seq[i+1]}{i+2}: moderate deamidation risk",
            })

    return sites
```

### 2.2 Oxidation Sites (Met, Trp, Free Cys)

```python
def find_oxidation_sites(sequence):
    """Find oxidation-prone residues.

    Met: sulfoxide formation (most common)
    Trp: various oxidation products
    Cys: disulfide scrambling if unpaired

    Returns:
        list of dicts with position, residue, risk
    """
    sites = []
    seq = sequence.upper()

    for i, aa in enumerate(seq):
        if aa == "M":
            # Met oxidation — always flag
            sites.append({
                "position": i + 1,
                "residue": "Met",
                "risk": "MODERATE",
                "note": f"Met{i+1}: sulfoxide risk (common in storage/process)",
            })
        elif aa == "W":
            # Trp oxidation — flag if surface-exposed (can't tell from sequence alone)
            sites.append({
                "position": i + 1,
                "residue": "Trp",
                "risk": "LOW",
                "note": f"Trp{i+1}: oxidation possible if surface-exposed",
            })
        elif aa == "C":
            sites.append({
                "position": i + 1,
                "residue": "Cys",
                "risk": "HIGH" if seq.count("C") % 2 != 0 else "MODERATE",
                "note": f"Cys{i+1}: check disulfide pairing",
            })

    return sites
```

### 2.3 Isomerization (Asp-X)

```python
def find_isomerization_sites(sequence):
    """Find Asp isomerization-prone sites.

    Asp-Gly and Asp-Ser are prone to iso-Asp formation.

    Returns:
        list of dicts
    """
    _ISO_HIGH = {"DG", "DS", "DT", "DN"}
    sites = []
    seq = sequence.upper()

    for i in range(len(seq) - 1):
        if seq[i] != "D":
            continue
        dipeptide = seq[i:i+2]
        if dipeptide in _ISO_HIGH:
            sites.append({
                "position": i + 1,
                "motif": dipeptide,
                "risk": "MODERATE",
                "note": f"Asp{i+1}-{seq[i+1]}{i+2}: iso-Asp formation risk",
            })

    return sites
```

### 2.4 Glycosylation Motifs (N-X-S/T)

```python
def find_glycosylation_sites(sequence):
    """Find N-linked glycosylation sequons (N-X-S/T, X≠P).

    Relevant for mammalian expression (CHO, HEK).
    Not relevant for E. coli (no glycosylation machinery).

    Returns:
        list of dicts
    """
    sites = []
    seq = sequence.upper()

    for i in range(len(seq) - 2):
        if seq[i] != "N":
            continue
        if seq[i+1] == "P":
            continue  # Pro blocks glycosylation
        if seq[i+2] in ("S", "T"):
            sites.append({
                "position": i + 1,
                "motif": seq[i:i+3],
                "risk": "INFO",
                "note": (f"N{i+1}-{seq[i+1]}{i+2}-{seq[i+2]}{i+3}: "
                        f"N-glycosylation sequon (relevant for mammalian expression)"),
            })

    return sites
```

## Phase 3: Aggregation & Solubility Prediction

### 3.1 Hydrophobic Patch Detection

```python
#!/opt/conda/envs/prot/bin/python
"""Detect hydrophobic patches in protein sequences."""
import numpy as np

# Kyte-Doolittle hydrophobicity scale
_KD_SCALE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "D": -3.5,
    "E": -3.5, "N": -3.5, "Q": -3.5, "K": -3.9, "R": -4.5,
}

def find_hydrophobic_patches(sequence, window=7, threshold=1.5, min_patch_len=5):
    """Find regions with sustained high hydrophobicity.

    Uses sliding-window Kyte-Doolittle average.
    Patches of sustained hydrophobicity (>threshold for >min_patch_len)
    indicate aggregation-prone regions.

    Args:
        sequence: amino acid string
        window: sliding window size
        threshold: KD average above this = hydrophobic
        min_patch_len: minimum consecutive positions above threshold

    Returns:
        list of patches with start, end, mean_hydrophobicity
    """
    seq = sequence.upper()
    L = len(seq)

    if L < window:
        return []

    # Compute sliding window average
    kd_values = np.array([_KD_SCALE.get(aa, 0) for aa in seq])
    kd_smooth = np.convolve(kd_values, np.ones(window)/window, mode="valid")

    # Find patches above threshold
    above = kd_smooth > threshold
    patches = []
    start = None

    for i, val in enumerate(above):
        if val and start is None:
            start = i
        elif not val and start is not None:
            if i - start >= min_patch_len:
                patches.append({
                    "start": start + 1,  # 1-indexed
                    "end": start + window + (i - start) - 1,
                    "length": i - start,
                    "mean_kd": float(kd_smooth[start:i].mean()),
                    "risk": "HIGH" if kd_smooth[start:i].mean() > 2.0 else "MODERATE",
                })
            start = None

    # Handle patch at end
    if start is not None and len(above) - start >= min_patch_len:
        patches.append({
            "start": start + 1,
            "end": L,
            "length": len(above) - start,
            "mean_kd": float(kd_smooth[start:].mean()),
            "risk": "HIGH" if kd_smooth[start:].mean() > 2.0 else "MODERATE",
        })

    return patches
```

### 3.2 Solubility Estimation

```python
def solubility_estimate(sequence):
    """Quick solubility risk assessment from sequence features.

    Based on empirical correlations:
    - Net charge magnitude: higher = more soluble
    - Fraction charged residues: higher = more soluble
    - GRAVY: more negative = more soluble
    - Hydrophobic patches: more/bigger = less soluble

    Returns:
        dict with solubility risk assessment
    """
    props = sequence_properties(sequence)
    patches = find_hydrophobic_patches(sequence)

    score = 0  # higher = better solubility
    reasons = []

    # Charge
    abs_charge = abs(props["charge_pH7"])
    if abs_charge > 5:
        score += 2
        reasons.append(f"Good net charge ({props['charge_pH7']:.1f})")
    elif abs_charge > 2:
        score += 1
    else:
        score -= 1
        reasons.append(f"Low net charge ({props['charge_pH7']:.1f}) — aggregation risk")

    # Fraction charged
    if props["frac_charged"] > 0.25:
        score += 1
        reasons.append(f"Good charged fraction ({props['frac_charged']:.0%})")
    elif props["frac_charged"] < 0.15:
        score -= 1
        reasons.append(f"Low charged fraction ({props['frac_charged']:.0%})")

    # GRAVY
    if props["gravy"] < -0.5:
        score += 1
        reasons.append("Hydrophilic (GRAVY < -0.5)")
    elif props["gravy"] > 0:
        score -= 2
        reasons.append(f"Hydrophobic (GRAVY={props['gravy']:.2f}) — high aggregation risk")

    # Hydrophobic patches
    n_patches = len(patches)
    if n_patches == 0:
        score += 1
        reasons.append("No hydrophobic patches")
    elif n_patches <= 2:
        reasons.append(f"{n_patches} hydrophobic patch(es) — moderate risk")
    else:
        score -= 2
        reasons.append(f"{n_patches} hydrophobic patches — high aggregation risk")

    # Classification
    if score >= 3:
        risk = "LOW"
    elif score >= 0:
        risk = "MODERATE"
    else:
        risk = "HIGH"

    return {
        "solubility_score": score,
        "solubility_risk": risk,
        "reasons": reasons,
        "n_hydrophobic_patches": n_patches,
    }
```

## Phase 4: Expression System Compatibility

### 4.1 E. coli Compatibility

```python
def ecoli_compatibility(sequence):
    """Check E. coli expression compatibility.

    Flags:
    - Free cysteines (reducing cytoplasm disrupts disulfides)
    - Rare codons (not checked here — need DNA sequence)
    - Large size (>60 kDa harder)
    - Glycosylation sites (won't be glycosylated)
    - Transmembrane regions (need special strains)

    Returns:
        dict with flags and recommendations
    """
    seq = sequence.upper()
    props = sequence_properties(seq)
    flags = []
    recs = []

    # Size
    if props["mw_da"] > 60000:
        flags.append(f"LARGE_PROTEIN ({props['mw_da']/1000:.0f} kDa) — may need optimization")
        recs.append("Consider: split into domains, or use pET-SUMO fusion")
    elif props["mw_da"] < 10000:
        flags.append(f"SMALL_PROTEIN ({props['mw_da']/1000:.1f} kDa) — may be unstable/degraded")
        recs.append("Consider: MBP or GST fusion for stabilization")

    # Cysteines
    n_cys = seq.count("C")
    if n_cys > 0:
        if n_cys % 2 == 0:
            flags.append(f"DISULFIDE_BONDS_LIKELY ({n_cys} Cys) — need oxidizing environment")
            recs.append("Use SHuffle or Origami strains, or secrete to periplasm")
        else:
            flags.append(f"ODD_CYSTEINES ({n_cys}) — unpaired Cys may cause aggregation")
            recs.append("Consider: mutate unpaired Cys to Ser, or use reducing conditions")

    # Glycosylation sites (won't work in E. coli)
    glyco_sites = find_glycosylation_sites(seq)
    if glyco_sites:
        flags.append(f"GLYCOSYLATION_SITES ({len(glyco_sites)}) — not glycosylated in E. coli")
        recs.append("If glycosylation needed, switch to CHO/HEK; otherwise OK to ignore")

    # Instability
    if props["instability_index"] > 40:
        flags.append("PREDICTED_UNSTABLE — may degrade rapidly")
        recs.append("Consider: lower induction temperature (16-18°C), or MBP fusion")

    return {
        "system": "E. coli",
        "flags": flags,
        "recommendations": recs,
        "n_flags": len(flags),
    }
```

### 4.2 Mammalian (CHO/HEK) Compatibility

```python
def mammalian_compatibility(sequence):
    """Check mammalian expression compatibility.

    Flags:
    - Free cysteines (oxidizing ER forms unwanted disulfides)
    - Missing signal peptide (for secreted proteins)
    - Aggregation-prone regions
    - N-glycosylation sites (will be glycosylated — may affect function)

    Returns:
        dict with flags and recommendations
    """
    seq = sequence.upper()
    props = sequence_properties(seq)
    flags = []
    recs = []

    # Cysteines
    n_cys = seq.count("C")
    if n_cys % 2 != 0:
        flags.append(f"ODD_CYSTEINES ({n_cys}) — unwanted disulfide formation in ER")
        recs.append("Mutate unpaired Cys to Ser before mammalian expression")

    # Glycosylation (will happen in CHO/HEK)
    glyco_sites = find_glycosylation_sites(seq)
    if glyco_sites:
        flags.append(f"GLYCOSYLATION_SITES ({len(glyco_sites)}) — will be glycosylated")
        recs.append("Verify: glycosylation at these sites acceptable for your application?")

    # Aggregation
    sol = solubility_estimate(seq)
    if sol["solubility_risk"] == "HIGH":
        flags.append("HIGH_AGGREGATION_RISK — may form aggregates in secretory pathway")
        recs.append("Consider: add solubility-enhancing mutations (charged surface residues)")

    return {
        "system": "mammalian (CHO/HEK)",
        "flags": flags,
        "recommendations": recs,
        "n_flags": len(flags),
    }
```

## Phase 5: Full Developability Report

### 5.1 Generate Red-Light Table

```python
#!/opt/conda/envs/prot/bin/python
"""Generate comprehensive developability report."""
import pandas as pd

def developability_report(sequence, name="design",
                           expression_system="ecoli"):
    """Run all developability checks and produce a red-light table.

    Args:
        sequence: amino acid string
        name: design name/label
        expression_system: "ecoli" or "mammalian"

    Returns:
        dict with all results + red_lights DataFrame
    """
    seq = sequence.upper().replace("X", "").replace("-", "")

    # Phase 1: Sequence properties
    props = sequence_properties(seq)
    comp_flags = composition_flags(seq)

    # Phase 2: Chemical liabilities
    deamidation = find_deamidation_sites(seq)
    oxidation = find_oxidation_sites(seq)
    isomerization = find_isomerization_sites(seq)
    glycosylation = find_glycosylation_sites(seq)

    # Phase 3: Aggregation
    patches = find_hydrophobic_patches(seq)
    solubility = solubility_estimate(seq)

    # Phase 4: Expression compatibility
    if expression_system == "ecoli":
        expr_compat = ecoli_compatibility(seq)
    else:
        expr_compat = mammalian_compatibility(seq)

    # Build red-light table
    red_lights = []

    # From properties
    for flag in props.get("flags", []):
        red_lights.append({"category": "biophysics", "issue": flag, "severity": "WARNING"})
    for flag in comp_flags:
        sev = "RED" if "ODD_CYSTEINE" in flag or "HOMOPOLYMER" in flag else "WARNING"
        red_lights.append({"category": "composition", "issue": flag, "severity": sev})

    # From chemical liabilities
    for site in deamidation:
        sev = "RED" if site["risk"] == "HIGH" else "WARNING"
        red_lights.append({"category": "deamidation", "issue": site["note"], "severity": sev})
    for site in oxidation:
        sev = "RED" if site["risk"] == "HIGH" else "WARNING"
        red_lights.append({"category": "oxidation", "issue": site["note"], "severity": sev})
    for site in isomerization:
        red_lights.append({"category": "isomerization", "issue": site["note"], "severity": "WARNING"})
    for site in glycosylation:
        red_lights.append({"category": "glycosylation", "issue": site["note"], "severity": "INFO"})

    # From aggregation
    for patch in patches:
        sev = "RED" if patch["risk"] == "HIGH" else "WARNING"
        red_lights.append({
            "category": "aggregation",
            "issue": f"Hydrophobic patch {patch['start']}-{patch['end']} "
                     f"(mean KD={patch['mean_kd']:.1f})",
            "severity": sev,
        })
    if solubility["solubility_risk"] == "HIGH":
        red_lights.append({
            "category": "solubility",
            "issue": "HIGH aggregation/solubility risk",
            "severity": "RED",
        })

    # From expression compatibility
    for flag in expr_compat.get("flags", []):
        sev = "RED" if "ODD_CYSTEINE" in flag else "WARNING"
        red_lights.append({"category": "expression", "issue": flag, "severity": sev})

    red_lights_df = pd.DataFrame(red_lights)

    # Summary
    n_red = len([r for r in red_lights if r["severity"] == "RED"])
    n_warn = len([r for r in red_lights if r["severity"] == "WARNING"])
    n_info = len([r for r in red_lights if r["severity"] == "INFO"])

    result = {
        "name": name,
        "sequence": seq,
        "length": len(seq),
        "properties": props,
        "deamidation_sites": deamidation,
        "oxidation_sites": oxidation,
        "isomerization_sites": isomerization,
        "glycosylation_sites": glycosylation,
        "hydrophobic_patches": patches,
        "solubility": solubility,
        "expression_compat": expr_compat,
        "red_lights": red_lights_df,
        "n_red": n_red,
        "n_warning": n_warn,
        "n_info": n_info,
        "expression_system": expression_system,
    }

    print(f"\n{'='*50}")
    print(f"Developability: {name} ({len(seq)} aa, {expression_system})")
    print(f"  RED flags:     {n_red}")
    print(f"  Warnings:      {n_warn}")
    print(f"  Info:          {n_info}")
    if n_red > 0:
        print(f"  ⛔ ACTION REQUIRED — {n_red} red flag(s) must be addressed")
    elif n_warn > 2:
        print(f"  ⚠️  Multiple warnings — review before proceeding")
    else:
        print(f"  ✅ Acceptable for experimental validation")
    print(f"{'='*50}")

    return result
```

### 5.2 Batch Developability Screen

```python
def batch_developability(sequences, names, expression_system="ecoli"):
    """Screen multiple designs for developability.

    Args:
        sequences: list of amino acid strings
        names: list of labels
        expression_system: "ecoli" or "mammalian"

    Returns:
        summary DataFrame, list of full reports
    """
    reports = []
    summaries = []

    for seq, name in zip(sequences, names):
        report = developability_report(seq, name, expression_system)
        reports.append(report)
        summaries.append({
            "name": name,
            "length": report["length"],
            "n_red": report["n_red"],
            "n_warning": report["n_warning"],
            "pI": report["properties"]["pI"],
            "gravy": report["properties"]["gravy"],
            "solubility_risk": report["solubility"]["solubility_risk"],
            "n_deamidation_high": sum(1 for s in report["deamidation_sites"] if s["risk"] == "HIGH"),
            "n_hydrophobic_patches": len(report["hydrophobic_patches"]),
            "n_free_cys": report["properties"].get("n_free_cys", 0),
        })

    summary_df = pd.DataFrame(summaries)
    summary_df = summary_df.sort_values("n_red")  # fewest red flags first

    return summary_df, reports
```

## Phase 6: Hard Gates

### 6.1 Absolute Hard Gates

```python
def developability_hard_gates(report):
    """Apply hard gates to a developability report.

    RED flags are hard fails. Designs with RED flags should not
    proceed to experiments without addressing the flagged issues.

    Args:
        report: from developability_report()

    Returns:
        pass/fail, list of blocking issues, list of fixable issues
    """
    blocking = []
    fixable = []

    for _, row in report["red_lights"].iterrows():
        if row["severity"] == "RED":
            # Check if fixable
            if "deamidation" in row["category"]:
                fixable.append(f"{row['issue']} → fix: mutate Asn to Gln")
            elif "ODD_CYSTEINE" in row["issue"]:
                fixable.append(f"{row['issue']} → fix: mutate Cys to Ser")
            elif "HOMOPOLYMER" in row["issue"]:
                blocking.append(row["issue"])
            else:
                blocking.append(row["issue"])

    passed = len(blocking) == 0

    return {
        "passed": passed,
        "n_blocking": len(blocking),
        "n_fixable": len(fixable),
        "blocking_issues": blocking,
        "fixable_issues": fixable,
        "recommendation": (
            "PASS — proceed to experiments" if passed and not fixable
            else "PASS_WITH_FIXES — address fixable issues first" if passed
            else "FAIL — blocking issues must be resolved"
        ),
    }
```

### 6.2 WT-Relative Mode (Compare Design vs Wildtype)

```python
#!/opt/conda/envs/prot/bin/python
"""WT-relative developability: only flag issues the design ADDS beyond WT.

Why this exists:
- IL-2 (1Z92) has ODD_CYSTEINES(3), instability>40, deamidation motifs
- These are inherent to the protein, not introduced by the design
- Absolute gates would reject EVERY design including the WT itself
- WT-relative mode: only veto if design is WORSE than WT

When to use:
- Always preferred when you have the WT sequence
- Essential for proteins with known developability quirks
  (disulfide-rich, intrinsically disordered regions, etc.)
"""

def developability_hard_gates_relative(design_report, wt_report):
    """Compare design developability against wildtype baseline.

    Only flags that are NEW in the design (not present in WT) are
    treated as blocking. Flags shared with WT are downgraded to INFO.

    Args:
        design_report: from developability_report(design_seq)
        wt_report: from developability_report(wt_seq)

    Returns:
        dict with relative assessment
    """
    # Extract flag sets
    wt_flags = set()
    if len(wt_report["red_lights"]) > 0:
        for _, row in wt_report["red_lights"].iterrows():
            wt_flags.add((row["category"], row["issue"]))

    design_flags = set()
    if len(design_report["red_lights"]) > 0:
        for _, row in design_report["red_lights"].iterrows():
            design_flags.add((row["category"], row["issue"]))

    # New flags: in design but NOT in WT
    new_flags = design_flags - wt_flags
    # Shared flags: in both design and WT
    shared_flags = design_flags & wt_flags
    # Resolved flags: in WT but NOT in design (design is better!)
    resolved_flags = wt_flags - design_flags

    # Only new RED flags are blocking
    blocking = []
    fixable = []
    for cat, issue in new_flags:
        # Check severity from design report
        match = design_report["red_lights"][
            (design_report["red_lights"]["category"] == cat) &
            (design_report["red_lights"]["issue"] == issue)
        ]
        if len(match) > 0 and match.iloc[0]["severity"] == "RED":
            if "deamidation" in cat:
                fixable.append(f"NEW: {issue} → fix: Asn→Gln")
            elif "ODD_CYSTEINE" in issue:
                fixable.append(f"NEW: {issue} → fix: Cys→Ser")
            else:
                blocking.append(f"NEW: {issue}")

    passed = len(blocking) == 0

    result = {
        "mode": "wt_relative",
        "passed": passed,
        "n_new_flags": len(new_flags),
        "n_shared_flags": len(shared_flags),
        "n_resolved_flags": len(resolved_flags),
        "n_blocking": len(blocking),
        "n_fixable": len(fixable),
        "blocking_issues": blocking,
        "fixable_issues": fixable,
        "shared_with_wt": [f"{cat}: {issue}" for cat, issue in shared_flags],
        "resolved_vs_wt": [f"{cat}: {issue}" for cat, issue in resolved_flags],
    }

    # Summary
    print(f"\n{'='*50}")
    print(f"Developability (WT-relative mode)")
    print(f"  Shared with WT (OK):  {len(shared_flags)}")
    print(f"  New in design:        {len(new_flags)}")
    print(f"  Resolved vs WT:       {len(resolved_flags)}")
    if blocking:
        print(f"  ⛔ NEW blocking:      {len(blocking)}")
    elif fixable:
        print(f"  ⚠️  NEW fixable:      {len(fixable)}")
    else:
        print(f"  ✅ No new issues vs WT")
    print(f"{'='*50}")

    if not blocking and not fixable:
        result["recommendation"] = ("PASS — design has no new developability "
                                    "issues beyond wildtype")
    elif not blocking:
        result["recommendation"] = ("PASS_WITH_FIXES — design has new fixable "
                                    "issues not in wildtype")
    else:
        result["recommendation"] = ("FAIL — design introduces new blocking "
                                    "issues not present in wildtype")

    return result
```

### 6.3 Quantitative WT Comparison

```python
def compare_properties_vs_wt(design_report, wt_report):
    """Compare biophysical properties: design vs WT.

    Flags significant deviations.

    Returns:
        dict with property deltas and flags
    """
    d = design_report["properties"]
    w = wt_report["properties"]

    deltas = {
        "delta_pI": d["pI"] - w["pI"],
        "delta_charge_pH7": d["charge_pH7"] - w["charge_pH7"],
        "delta_gravy": d["gravy"] - w["gravy"],
        "delta_instability": d["instability_index"] - w["instability_index"],
    }

    flags = []
    if abs(deltas["delta_pI"]) > 1.0:
        flags.append(f"pI shifted by {deltas['delta_pI']:+.1f} "
                     f"({w['pI']:.1f}→{d['pI']:.1f})")
    if abs(deltas["delta_charge_pH7"]) > 3:
        flags.append(f"Net charge shifted by {deltas['delta_charge_pH7']:+.1f}")
    if deltas["delta_gravy"] > 0.3:
        flags.append(f"Design more hydrophobic (ΔGRAVY={deltas['delta_gravy']:+.2f})")
    if deltas["delta_instability"] > 10:
        flags.append(f"Instability increased (Δ={deltas['delta_instability']:+.1f})")

    if not flags:
        flags.append("Properties similar to WT — no concerns")

    deltas["flags"] = flags
    return deltas
```

### 6.4 When to Use Which Mode

```
Which developability gate mode?
│
├── Do you have the WT sequence?
│   ├── YES → Use WT-relative mode (Phase 6.2)
│   │   └── Also run quantitative comparison (Phase 6.3)
│   └── NO → Use absolute mode (Phase 6.1)
│
└── Design has same flags as WT?
    ├── YES → Shared flags are OK (protein's inherent properties)
    │         Only NEW flags matter
    └── NO (design has new flags) → Investigate and fix/reject
```

## Phase 7: Integration with prot-dbtl

```
prot-dbtl pipeline position:
│
Step 1: prot-structure (quality check)
Step 2: prot-interface (hotspots + constraints)
Step 3: prot-seqdesign (design)
Step 4: prot-stability (LLR + consensus)
Step 5: prot-structure (fold-back)
Step 6: prot-interface (contact preservation)
│
▼
Step 7: prot-developability (RED-LIGHT CHECK)  ← NEW
│
├── All candidates pass? → Proceed to DNA ordering
├── Some have fixable issues? → Apply fixes, re-run fold-back
└── Blocking issues? → Reject or redesign with constraints
```

## Failure Modes

1. **Skipping developability screening.** The most common failure. "It passed all computational gates" means nothing if it has 3 unpaired cysteines and a hydrophobic patch. Takes 1 second to check.

2. **Treating all red flags equally.** An unpaired cysteine is fixable (C→S). A massive hydrophobic core is not. Distinguish blocking vs fixable issues.

3. **Ignoring expression system.** A free Cys is fine in reducing E. coli cytoplasm but catastrophic in CHO secretory pathway. Always specify the system.

4. **Over-engineering for developability.** Don't mutate every deamidation site — some are buried and never deamidate. Use structural context (prot-structure) to distinguish exposed vs buried liabilities.

5. **Not checking the wildtype.** If the wildtype has 5 deamidation sites and expresses fine, don't penalize designs for having 5 deamidation sites. Compare to wildtype.

## One-Sentence Rule

**Run every candidate through the red-light table before ordering DNA, always compare against wildtype baseline (only NEW flags matter), and remember that free cysteines, deamidation hotspots, and aggregation patches kill more designs than bad scores ever will.**
