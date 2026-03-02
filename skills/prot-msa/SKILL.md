---
name: prot-msa
description: Multiple sequence alignment, conservation analysis, and position-specific mutation risk scoring. Uses mmseqs2/hhblits for real MSA when available, with ESM-2 pseudo-conservation as zero-dependency fallback. Provides family-calibrated conservation scores, consensus sequences, coevolution signals, and low-risk mutation suggestions. Triggers whenever you need to assess how conserved a position is, build a consensus sequence, identify coevolving positions, or calibrate ESM-2 LLR with family context.
homepage: https://github.com/soedinglab/MMseqs2
metadata: { "openclaw": { "emoji": "📊", "requires": { "bins": ["python3"], "python": ["numpy", "pandas", "biopython", "torch", "fair-esm"] } } }
---

# Protein MSA: Conservation, Consensus & Family-Calibrated Design

ESM-2 LLR is a broad evolutionary prior trained on all proteins. This skill adds **family-specific** context: how conserved is each position in YOUR protein's family? A mutation that ESM-2 calls "unnatural" might be perfectly normal in your specific family — and vice versa.

## When to Use

- **Calibrating ESM-2 LLR**: is a low LLR score because the mutation is truly bad, or because ESM-2's global prior is too strict for this family?
- **Building a consensus sequence** for stability engineering
- **Identifying coevolving positions** that should be mutated together
- **Risk-stratifying mutations**: family-conserved positions are higher risk than variable ones
- **Generating "safe mutation" lists**: positions where the family already shows variation
- **Informing prot-dbtl constraints**: freeze conserved positions, redesign variable ones

## Core Philosophy

1. **Family context beats global prior.** ESM-2 sees all proteins; your protein lives in a specific family. A position that's 100% conserved in the family is truly important, even if ESM-2 thinks a mutation there is OK. Conversely, a position that varies across the family is safe to mutate, even if ESM-2 penalizes it.
2. **Conservation is position-specific, not sequence-level.** Don't say "this protein is conserved." Say "position 45 is 98% conserved (always Trp), position 72 is variable (6 different amino acids in the family)."
3. **Consensus mutations are low-risk stabilizers.** If position 45 is Ala in your sequence but Ser in 80% of homologs, mutating to Ser is likely stabilizing. This is the basis of consensus design.
4. **Coevolution reveals structural constraints.** If positions 30 and 65 always change together (coevolve), mutating one without the other may destabilize the protein.
5. **When you can't build a real MSA, ESM-2 provides a reasonable proxy.** ESM-2's per-position probability distribution is essentially a "virtual MSA" learned from millions of sequences. Not as specific as a real MSA, but much better than nothing.

## Phase 0: Environment Verification

```python
#!/opt/conda/envs/prot/bin/python
"""Phase 0: Check MSA tools availability."""
import importlib, shutil, subprocess

print("=== prot-msa Phase 0 ===\n")

# Core dependencies
deps = {"numpy": None, "pandas": None, "Bio": "biopython", "torch": None, "esm": "fair-esm"}
missing = []
for mod, pkg in deps.items():
    try:
        importlib.import_module(mod)
    except ImportError:
        missing.append(pkg or mod)

if missing:
    print(f"❌ Missing: {missing}")
else:
    print("✅ Core dependencies OK")

# MSA tools (optional but preferred)
msa_backend = "ESM2_PSEUDO"  # fallback default

# Check mmseqs2
if shutil.which("mmseqs"):
    print("✅ mmseqs2 found — real MSA available")
    msa_backend = "MMSEQS2"
else:
    print("ℹ️  mmseqs2 not found — will use ESM-2 pseudo-conservation")
    print("   Install: conda install -c bioconda mmseqs2")

# Check hhblits (alternative)
if shutil.which("hhblits"):
    print("✅ hhblits found — HHblits MSA available")
    if msa_backend == "ESM2_PSEUDO":
        msa_backend = "HHBLITS"
else:
    print("ℹ️  hhblits not found")

# Check databases
import os
UNIREF_DB = os.environ.get("UNIREF_DB", "")
if UNIREF_DB and os.path.exists(UNIREF_DB):
    print(f"✅ UniRef database: {UNIREF_DB}")
else:
    print("ℹ️  No UniRef database set (set UNIREF_DB env var)")
    if msa_backend in ("MMSEQS2", "HHBLITS"):
        print("   ⚠️  MSA tools found but no database — falling back to ESM-2")
        msa_backend = "ESM2_PSEUDO"

# GPU
import torch
if torch.cuda.is_available():
    vram = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f"✅ GPU: {torch.cuda.get_device_name(0)} ({vram:.1f} GB)")

print(f"\n→ MSA_BACKEND = {msa_backend}")
print("\n=== Phase 0 complete ===")
```

## Phase 1: Build MSA (When Tools Available)

### 1.1 mmseqs2 MSA

```bash
# Build MSA with mmseqs2 (fast, good for large databases)
# Requires: mmseqs2 installed + UniRef30/UniRef50 database

mmseqs easy-search \
  query.fasta \
  $UNIREF_DB \
  result.m8 \
  tmp/ \
  --max-seqs 1000 \
  -e 1e-5

# Convert to FASTA alignment
mmseqs easy-msa \
  query.fasta \
  $UNIREF_DB \
  msa_output \
  tmp/ \
  --max-seqs 500
```

### 1.2 Parse MSA and Compute Conservation

```python
#!/opt/conda/envs/prot/bin/python
"""Parse MSA and compute per-position conservation."""
import numpy as np
import pandas as pd
from Bio import AlignIO
from collections import Counter

def parse_msa(msa_path, format="fasta"):
    """Parse a multiple sequence alignment file.

    Args:
        msa_path: path to MSA file (FASTA, Stockholm, etc.)
        format: Biopython alignment format

    Returns:
        alignment object, list of sequences, query sequence
    """
    alignment = AlignIO.read(msa_path, format)
    sequences = [str(record.seq) for record in alignment]
    query = sequences[0]  # first sequence is usually the query

    print(f"MSA: {len(sequences)} sequences, {len(query)} positions")
    return alignment, sequences, query


def compute_conservation(sequences, query=None):
    """Compute per-position conservation from MSA.

    Metrics per position:
    - frequency: fraction of most common amino acid
    - entropy: Shannon entropy (0 = fully conserved, ~4.3 = uniform)
    - n_types: number of distinct amino acids observed
    - consensus_aa: most common amino acid
    - gap_fraction: fraction of gaps at this position

    Args:
        sequences: list of aligned sequences (same length)
        query: query sequence (first in list if None)

    Returns:
        DataFrame with per-position conservation metrics
    """
    if query is None:
        query = sequences[0]

    L = len(query)
    n_seqs = len(sequences)

    rows = []
    for i in range(L):
        col = [s[i] for s in sequences if i < len(s)]
        col_no_gap = [aa for aa in col if aa not in ("-", ".", "X")]

        if not col_no_gap:
            rows.append({
                "position": i + 1, "query_aa": query[i],
                "consensus_aa": "-", "frequency": 0, "entropy": 0,
                "n_types": 0, "gap_fraction": 1.0,
            })
            continue

        counts = Counter(col_no_gap)
        total = sum(counts.values())

        # Most common
        consensus_aa, max_count = counts.most_common(1)[0]
        frequency = max_count / total

        # Shannon entropy
        probs = np.array([c / total for c in counts.values()])
        entropy = -np.sum(probs * np.log2(probs + 1e-10))

        # Gap fraction
        n_gaps = sum(1 for aa in col if aa in ("-", "."))
        gap_fraction = n_gaps / len(col)

        rows.append({
            "position": i + 1,
            "query_aa": query[i],
            "consensus_aa": consensus_aa,
            "frequency": float(frequency),
            "entropy": float(entropy),
            "n_types": len(counts),
            "gap_fraction": float(gap_fraction),
        })

    df = pd.DataFrame(rows)

    # Conservation classification
    df["conservation"] = df["entropy"].apply(_classify_conservation)

    return df


def _classify_conservation(entropy):
    """Classify conservation level from Shannon entropy.

    Lower entropy = more conserved.
    """
    if entropy < 0.5:
        return "HIGHLY_CONSERVED"   # almost invariant
    elif entropy < 1.5:
        return "CONSERVED"          # strong preference
    elif entropy < 2.5:
        return "MODERATE"           # some variation
    else:
        return "VARIABLE"           # many amino acids observed
```

## Phase 2: ESM-2 Pseudo-Conservation (Zero-Dependency Fallback)

### 2.1 ESM-2 Per-Position Entropy

```python
#!/opt/conda/envs/prot/bin/python
"""Compute pseudo-conservation from ESM-2 per-position distributions.

ESM-2's per-position probability distribution over amino acids is
essentially a "virtual MSA" learned from millions of sequences.
Lower entropy = ESM-2 thinks this position is more constrained.
"""
import torch
import numpy as np
import pandas as pd

def esm2_pseudo_conservation(sequence, model=None, alphabet=None,
                              batch_converter=None, device=None):
    """Compute pseudo-conservation scores from ESM-2.

    For each position, computes:
    - Predicted probability distribution over 20 standard amino acids
    - Shannon entropy of that distribution
    - Rank of the query amino acid
    - Top-3 predicted amino acids (consensus proxy)

    Args:
        sequence: amino acid string

    Returns:
        DataFrame with per-position pseudo-conservation metrics
    """
    if model is None:
        import esm
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        model = model.eval().to(device)
        batch_converter = alphabet.get_batch_converter()

    data = [("protein", sequence)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    with torch.no_grad():
        results = model(tokens, repr_layers=[], return_contacts=False)
        logits = results["logits"][0]  # (L+2, vocab)

    standard_aas = "ACDEFGHIKLMNPQRSTVWY"
    aa_indices = {aa: alphabet.get_idx(aa) for aa in standard_aas}

    rows = []
    for i, aa in enumerate(sequence):
        pos_logits = logits[i + 1]  # +1 for BOS
        probs_all = torch.softmax(pos_logits, dim=-1)

        # Extract probabilities for standard AAs
        probs = {res: float(probs_all[idx]) for res, idx in aa_indices.items()}
        prob_array = np.array(list(probs.values()))

        # Normalize to sum=1 (exclude special tokens)
        prob_array = prob_array / prob_array.sum()

        # Shannon entropy
        entropy = -np.sum(prob_array * np.log2(prob_array + 1e-10))

        # Rank of query AA
        ranked = sorted(probs.items(), key=lambda x: -x[1])
        rank = next(j for j, (r, _) in enumerate(ranked) if r == aa) + 1

        # Top-3 consensus
        top3 = [r[0] for r in ranked[:3]]
        query_in_top3 = aa in top3

        # Query probability
        query_prob = probs.get(aa, 0)

        rows.append({
            "position": i + 1,
            "query_aa": aa,
            "query_prob": float(query_prob),
            "query_rank": rank,
            "query_in_top3": query_in_top3,
            "entropy": float(entropy),
            "top3": "".join(top3),
            "consensus_aa": ranked[0][0],
            "consensus_prob": float(ranked[0][1]),
        })

    df = pd.DataFrame(rows)
    df["conservation"] = df["entropy"].apply(_classify_conservation)

    n_high = (df["conservation"] == "HIGHLY_CONSERVED").sum()
    n_var = (df["conservation"] == "VARIABLE").sum()
    print(f"ESM-2 pseudo-conservation: {n_high} highly conserved, "
          f"{n_var} variable out of {len(df)} positions")

    return df
```

### 2.2 Interpretation: ESM-2 Pseudo-Conservation vs Real MSA

| Metric | Real MSA | ESM-2 pseudo |
|--------|----------|-------------|
| Specificity | Family-specific | Global (all proteins) |
| Sensitivity to rare families | Excellent | Poor (may miss family-specific patterns) |
| Speed | Slow (needs database search) | Fast (single forward pass) |
| Dependencies | mmseqs2 + database | ESM-2 only |
| Use case | Final analysis, publication | Quick screening, no-database environments |

> **Rule of thumb:** If ESM-2 pseudo-conservation says VARIABLE but real MSA says CONSERVED → trust the MSA (family-specific signal). If ESM-2 says CONSERVED but MSA says VARIABLE → also trust the MSA (ESM-2 is being overly cautious from global patterns).

## Phase 3: Consensus Design

### 3.1 Build Consensus Sequence

```python
#!/opt/conda/envs/prot/bin/python
"""Build consensus sequence and identify stabilizing mutations."""

def consensus_design(query_sequence, conservation_df):
    """Identify positions where query differs from consensus.

    Consensus mutations (query → consensus AA) are generally stabilizing
    because they align with what evolution has selected for in this family.

    Args:
        query_sequence: your protein sequence
        conservation_df: from compute_conservation() or esm2_pseudo_conservation()

    Returns:
        DataFrame of suggested consensus mutations, sorted by likely impact
    """
    mutations = []

    for _, row in conservation_df.iterrows():
        pos = int(row["position"]) - 1  # 0-indexed
        query_aa = row["query_aa"]
        consensus_aa = row["consensus_aa"]

        if query_aa == consensus_aa:
            continue  # already at consensus

        if consensus_aa in ("-", ".", "X"):
            continue  # gap or unknown

        entropy = row["entropy"]
        cons = row["conservation"]

        # Only suggest mutations where there's clear consensus signal
        if cons in ("HIGHLY_CONSERVED", "CONSERVED"):
            risk = "LOW"
            rationale = f"Position strongly prefers {consensus_aa} (entropy={entropy:.2f})"
        elif cons == "MODERATE":
            risk = "MEDIUM"
            rationale = f"Position has moderate preference for {consensus_aa}"
        else:
            risk = "HIGH"
            rationale = f"Position is variable — consensus {consensus_aa} may not matter"

        mutations.append({
            "position": pos + 1,
            "query_aa": query_aa,
            "consensus_aa": consensus_aa,
            "mutation": f"{query_aa}{pos+1}{consensus_aa}",
            "entropy": float(entropy),
            "conservation": cons,
            "risk": risk,
            "rationale": rationale,
            "consensus_prob": float(row.get("consensus_prob", row.get("frequency", 0))),
        })

    df = pd.DataFrame(mutations)
    if len(df) > 0:
        df = df.sort_values("entropy")  # most conserved first = lowest risk

    n_low = (df["risk"] == "LOW").sum() if len(df) > 0 else 0
    n_med = (df["risk"] == "MEDIUM").sum() if len(df) > 0 else 0
    print(f"Consensus mutations: {n_low} low-risk, {n_med} medium-risk, "
          f"{len(df)-n_low-n_med} high-risk")

    return df
```

### 3.2 Consensus Mutations as prot-dbtl Constraints

```python
def consensus_to_constraints(consensus_mutations_df, chain_id, max_mutations=5):
    """Convert top consensus mutations into design constraints.

    Strategy: suggest a small number (3-5) of low-risk consensus mutations
    as a "stability boost" that can be combined with interface redesign.

    Args:
        consensus_mutations_df: from consensus_design()
        chain_id: chain ID
        max_mutations: max number to suggest

    Returns:
        dict with suggested mutations and reasoning
    """
    low_risk = consensus_mutations_df[consensus_mutations_df["risk"] == "LOW"]
    top = low_risk.head(max_mutations)

    suggestions = []
    for _, row in top.iterrows():
        suggestions.append({
            "position": f"{chain_id}{row['position']}",
            "mutation": row["mutation"],
            "rationale": row["rationale"],
        })

    print(f"\nTop {len(suggestions)} consensus mutations for stability boost:")
    for s in suggestions:
        print(f"  {s['mutation']}: {s['rationale']}")

    return {
        "suggested_mutations": suggestions,
        "strategy": "consensus_stabilization",
        "note": "Apply these alongside interface redesign to compensate for destabilizing interface mutations",
    }
```

## Phase 4: LLR Calibration with Conservation

### 4.1 Calibrate ESM-2 LLR Against Family Conservation

```python
#!/opt/conda/envs/prot/bin/python
"""Calibrate ESM-2 LLR using family conservation context.

This is the key integration with prot-stability:
- A mutation with LLR=-3 at a VARIABLE position → probably fine
- A mutation with LLR=-3 at a HIGHLY_CONSERVED position → probably bad
- A mutation with LLR=-1 at a HIGHLY_CONSERVED position → definitely bad
"""
import pandas as pd

def calibrate_llr_with_conservation(stability_df, conservation_df):
    """Add conservation context to stability LLR scores.

    Args:
        stability_df: per-mutation LLR from prot-stability
            Expected columns: position (1-indexed), wt_aa, llr
        conservation_df: from compute_conservation() or esm2_pseudo_conservation()
            Expected columns: position, conservation, entropy

    Returns:
        merged DataFrame with calibrated risk assessment
    """
    # Merge on position
    merged = stability_df.merge(
        conservation_df[["position", "conservation", "entropy", "consensus_aa",
                        "query_in_top3"]],
        left_on="position" if "position" in stability_df.columns else "position_1idx",
        right_on="position",
        how="left",
        suffixes=("", "_cons"),
    )

    # Calibrated risk
    risks = []
    for _, row in merged.iterrows():
        llr = row.get("llr", 0)
        cons = row.get("conservation", "MODERATE")

        if cons == "HIGHLY_CONSERVED":
            # Any negative LLR at a conserved position is concerning
            if llr < -0.5:
                risk = "HIGH"
                note = "Conserved position + unfavorable LLR → likely destabilizing"
            elif llr < 0:
                risk = "MODERATE"
                note = "Conserved position + mildly unfavorable → caution"
            else:
                risk = "LOW"
                note = "Conserved but LLR favorable → unusual but OK"
        elif cons == "CONSERVED":
            if llr < -2.0:
                risk = "HIGH"
                note = "Conserved position + strongly unfavorable LLR"
            elif llr < -0.5:
                risk = "MODERATE"
                note = "Conserved position + unfavorable → review"
            else:
                risk = "LOW"
                note = "Conserved and LLR acceptable"
        elif cons == "MODERATE":
            if llr < -3.0:
                risk = "HIGH"
                note = "Moderately conserved + very low LLR"
            elif llr < -1.0:
                risk = "MODERATE"
                note = "Some conservation + unfavorable LLR"
            else:
                risk = "LOW"
                note = "Moderate conservation, LLR acceptable"
        else:  # VARIABLE
            if llr < -4.0:
                risk = "MODERATE"
                note = "Variable position but extreme LLR — check for structural role"
            else:
                risk = "LOW"
                note = "Variable position — LLR probably too strict here"

        risks.append({"calibrated_risk": risk, "calibration_note": note})

    risk_df = pd.DataFrame(risks)
    merged = pd.concat([merged.reset_index(drop=True), risk_df], axis=1)

    # Summary
    if len(merged) > 0:
        n_high = (merged["calibrated_risk"] == "HIGH").sum()
        n_mod = (merged["calibrated_risk"] == "MODERATE").sum()
        n_low = (merged["calibrated_risk"] == "LOW").sum()
        print(f"Calibrated risk: {n_high} HIGH, {n_mod} MODERATE, {n_low} LOW")

    return merged
```

### 4.2 Calibrated Risk Thresholds

| Conservation | LLR < -0.5 | LLR -0.5 to -2 | LLR < -2 | LLR < -4 |
|-------------|------------|-----------------|----------|----------|
| **HIGHLY_CONSERVED** | HIGH | HIGH | HIGH | HIGH |
| **CONSERVED** | LOW | MODERATE | HIGH | HIGH |
| **MODERATE** | LOW | LOW | MODERATE | HIGH |
| **VARIABLE** | LOW | LOW | LOW | MODERATE |

> **Key insight:** At VARIABLE positions, even LLR=-3 is LOW risk because the family already tolerates variation there. At HIGHLY_CONSERVED positions, even LLR=-0.5 is HIGH risk because the position is under strong selection.

## Phase 5: Integration with Other Skills

### 5.1 Connection Map

```
prot-msa provides conservation context to:
│
├── prot-stability (Phase 4 calibration)
│   └── calibrate_llr_with_conservation() → more accurate risk assessment
│
├── prot-interface (constraint refinement)
│   └── Variable interface positions → safe to redesign
│   └── Conserved interface positions → likely structural, freeze
│
├── prot-seqdesign (informed constraints)
│   └── Consensus mutations → suggest as "stability boost"
│   └── Variable positions → higher temperature OK
│   └── Conserved positions → lower temperature or freeze
│
└── prot-dbtl (cycle planning)
    └── Conservation profile → initial constraint generation for Cycle 1
    └── Consensus mutations → combine with interface redesign for stability
```

### 5.2 Quick Integration Example

```python
#!/opt/conda/envs/prot/bin/python
"""Example: integrate conservation with stability + interface for 1BRS."""

def integrated_assessment(sequence, stability_per_mut_df, interface_burial_df,
                           msa_backend="ESM2_PSEUDO"):
    """Full assessment combining conservation + stability + burial.

    This is the comprehensive risk assessment that resolves the
    "fold-back PASS + interface PASS + LLR FAIL" conflict seen in 1BRS.

    Args:
        sequence: protein sequence
        stability_per_mut_df: from prot-stability score_design_vs_wildtype()
        interface_burial_df: from prot-interface annotate_interface_burial()
        msa_backend: "MMSEQS2", "HHBLITS", or "ESM2_PSEUDO"

    Returns:
        comprehensive per-mutation risk table
    """
    # Step 1: Get conservation
    if msa_backend == "ESM2_PSEUDO":
        cons_df = esm2_pseudo_conservation(sequence)
    else:
        # Would use real MSA here
        cons_df = esm2_pseudo_conservation(sequence)

    # Step 2: Calibrate LLR with conservation
    calibrated = calibrate_llr_with_conservation(stability_per_mut_df, cons_df)

    # Step 3: Add burial context (from prot-interface)
    if interface_burial_df is not None and len(interface_burial_df) > 0:
        burial_map = dict(zip(interface_burial_df["resseq"],
                             interface_burial_df["burial_rel"]))
        calibrated["burial"] = calibrated["position"].map(burial_map).fillna("UNKNOWN")
    else:
        calibrated["burial"] = "UNKNOWN"

    # Step 4: Final verdict
    verdicts = []
    for _, row in calibrated.iterrows():
        cal_risk = row.get("calibrated_risk", "MODERATE")
        burial = row.get("burial", "UNKNOWN")

        if cal_risk == "HIGH" and burial == "CORE":
            verdict = "REJECT — conserved + core + unfavorable LLR"
        elif cal_risk == "HIGH" and burial == "SURFACE":
            verdict = "REVIEW — conserved but surface, check if structural"
        elif cal_risk == "LOW" and burial in ("SURFACE", "PARTIAL"):
            verdict = "ACCEPT — variable/surface, LLR probably too strict"
        elif cal_risk == "MODERATE":
            verdict = "REVIEW — ambiguous, needs additional evidence"
        else:
            verdict = f"REVIEW — {cal_risk} risk, {burial} burial"

        verdicts.append(verdict)

    calibrated["verdict"] = verdicts

    return calibrated
```

## Phase 6: VRAM Budget & Timing

| Task | VRAM | Wall time | Notes |
|------|------|-----------|-------|
| ESM-2 pseudo-conservation | ~2 GB | ~2 sec | Single forward pass |
| mmseqs2 MSA (500 seqs) | 0 GB | ~30-60 sec | CPU, needs database |
| Conservation from MSA | 0 GB | ~1 sec | CPU, pure Python |
| Consensus design | 0 GB | ~1 sec | CPU |
| LLR calibration | 0 GB | ~1 sec | CPU, merge operation |

## Failure Modes

1. **Using ESM-2 pseudo-conservation as a substitute for real MSA.** ESM-2 captures global patterns but misses family-specific signals. If your protein is from a rare family or has unusual features, ESM-2 pseudo-conservation will be misleading. Always prefer real MSA when available.

2. **Over-trusting consensus.** Consensus mutations are LOW risk, not NO risk. The consensus amino acid might be preferred for a reason unrelated to stability (e.g., it's the ancestral state). Always check with fold-back validation.

3. **Ignoring gap-rich positions.** If >30% of MSA sequences have a gap at a position, the conservation signal is unreliable. Flag these positions as "uncertain" rather than "variable."

4. **Assuming all conserved positions are equal.** A position conserved as "always Trp" (aromatic core packing) is different from "always Asp" (catalytic residue). Conservation tells you the position matters; contact analysis (prot-interface) tells you WHY.

5. **Not recalculating after mutations.** Conservation is context-dependent. After introducing 5 mutations, the ESM-2 pseudo-conservation may change. For multi-mutation designs, recompute conservation in the mutant sequence context.

## One-Sentence Rule

**Use family conservation to calibrate ESM-2 LLR — variable positions tolerate mutations that LLR penalizes, conserved positions are dangerous even when LLR is mild — and suggest low-risk consensus mutations as stability anchors alongside functional redesign.**
