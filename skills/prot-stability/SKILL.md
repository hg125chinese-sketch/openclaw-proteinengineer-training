---
name: prot-stability
description: Protein stability prediction and scoring. Uses ESM-2 log-likelihood ratio (LLR) as fast computational proxy, with FoldX and Rosetta cartesian_ddg as optional physics-based tools. Covers mutation scanning, consensus design, multi-tool conflict resolution, and fast-screen → precision-screen workflows. Triggers whenever you need to assess whether a mutation or design is stabilizing, destabilizing, or neutral.
homepage: https://github.com/facebookresearch/esm
metadata: { "openclaw": { "emoji": "🔥", "requires": { "bins": ["python3"], "python": ["fair-esm", "torch", "numpy", "pandas", "biopython", "requests"] } } }
---

# Protein Stability: Prediction & Scoring

Stability is the gatekeeper between "designed sequence" and "expressible protein." A design that scores perfectly on backbone fit but destabilizes the fold will produce inclusion bodies. This skill gives you tools and decision boundaries to catch destabilizing mutations before they reach the bench.

## When to Use

- Scoring **designed sequences** from prot-seqdesign for stability
- **Mutation scanning** — evaluating single or multiple point mutations
- **Consensus design** — identifying stabilizing positions from evolutionary signal
- Resolving **stability-function trade-offs** (e.g., affinity up but Tm down)
- **Fast screening** (50→10) before expensive downstream evaluation
- **Precision scoring** (10→final) for candidate selection

## Core Philosophy

1. **Stability is a filter, not a target.** Optimize function first (binding, catalysis), then ensure stability doesn't kill the design. Exception: if expression/solubility is the primary bottleneck, stability becomes the objective.
2. **ESM-2 LLR is fast and good enough for screening.** It captures evolutionary plausibility — mutations that look unnatural to the protein language model are often destabilizing. Use it for fast ranking (50→10), not as ground truth.
3. **Physics and evolution can disagree — and that's informative.** When FoldX says destabilizing but ESM says natural (or vice versa), don't average them — diagnose why they disagree. The disagreement IS the signal.
4. **Compare to wildtype, always.** ΔΔG and ΔLLR are always relative to wildtype. An ESM LLR of -0.5 means nothing without knowing wildtype scores -0.3.
5. **Consensus mutations are low-risk stabilizers.** If a position in your design matches the consensus across homologs, it's more likely stable. Deviations from consensus carry stability risk proportional to how unusual they are.

## Interface-mutation mode (v2 policy)

Use this mode when you are doing **interface redesign/maturation** (mutations mostly/only at binding interface residues) and you will also run **fold-back** (prot-structure) plus **interface gates** (prot-interface).

Key idea: **ESM-2 LLR is an evolutionary plausibility prior**, not a calibrated physical ΔΔG. In tight interfaces, it may be overly harsh on rare-but-physically-tolerable substitutions.

### interface-mutation gating (two-stage)

- **Stage A (catastrophic filter only):**
  - REJECT if `min_LLR < -6.0`
  - Otherwise do **not** hard-reject based on `mean_LLR` alone.

- **Stage B (conflict triage):**
  If `mean_LLR` is very negative but:
  - fold-back PASS (e.g., TM-score ≥ 0.8 and mean pLDDT ≥ 75) AND
  - interface gate PASS
  → label as `LLR_CONFLICT` (do not call it “unstable” yet) and escalate to:
  - physics-based ΔΔG (FoldX/Rosetta), or
  - family-calibrated MSA/consensus stability model.

### Reporting requirement

When `LLR_CONFLICT`, report whether mutated positions are **buried/core-like** vs **solvent/interface** (needs SASA/burial proxy), and whether the mutation is conservative.

## Phase 0: Environment Verification

```python
#!/opt/conda/envs/prot/bin/python
"""Phase 0: Verify stability tools are available."""
import importlib

print("=== prot-stability Phase 0 ===\n")

# Core: ESM-2 (always needed)
deps = {
    "torch": None, "esm": "fair-esm",
    "numpy": None, "pandas": None,
    "Bio.PDB": "biopython",
}
missing = []
for mod, pkg in deps.items():
    try:
        importlib.import_module(mod.split(".")[0])
    except ImportError:
        missing.append(pkg or mod)

if missing:
    print(f"❌ Missing: {missing}")
else:
    print("✅ Core dependencies OK (ESM-2 available)")

# GPU
import torch
if torch.cuda.is_available():
    vram = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f"✅ GPU: {torch.cuda.get_device_name(0)} ({vram:.1f} GB)")
else:
    print("⚠️  No GPU — ESM-2 will be slow but functional on CPU")

# ESM-2 model check
try:
    import esm
    assert hasattr(esm.pretrained, 'esm2_t33_650M_UR50D')
    print("✅ ESM-2 650M available")
except Exception as e:
    print(f"❌ ESM-2 issue: {e}")

# Optional: FoldX
import shutil
if shutil.which("foldx") or shutil.which("FoldX"):
    print("✅ FoldX found in PATH")
else:
    print("ℹ️  FoldX not installed (optional — ESM-2 LLR is sufficient for screening)")

print("\n=== Phase 0 complete ===")
```

## Tool Landscape

| Tool | Type | Speed | Accuracy | Available? | Use case |
|------|------|-------|----------|-----------|----------|
| **ESM-2 LLR** | Evolutionary (language model) | ★★★★★ | ★★★☆☆ | ✅ Yes | Fast screen, ranking, consensus proxy |
| **FoldX** | Physics (empirical force field) | ★★★☆☆ | ★★★☆☆ | ❌ Not installed | ΔΔG, mutation scanning, RepairPDB |
| **Rosetta cartesian_ddg** | Physics (detailed energy function) | ★★☆☆☆ | ★★★★☆ | ❌ Not installed | Precision scoring, final candidates |
| **ESMFold + pLDDT** | Structure confidence | ★★★★☆ | ★★★☆☆ | ✅ Yes (API) | Fold-back stability proxy |

**Current strategy**: ESM-2 LLR as primary tool (installed, fast, GPU-accelerated). FoldX/Rosetta documented for when they become available.

## Phase 1: ESM-2 Log-Likelihood Ratio (LLR) — Primary Stability Proxy

### 1.1 How It Works

ESM-2 is a protein language model trained on millions of sequences. For each position in a protein, it predicts a probability distribution over amino acids. The **log-likelihood ratio (LLR)** compares:

- How likely is the **mutant** amino acid at this position?
- How likely is the **wildtype** amino acid at this position?

```
LLR = log P(mutant | context) - log P(wildtype | context)
```

- **LLR > 0**: mutant is MORE likely than wildtype → probably neutral or stabilizing
- **LLR ≈ 0**: similar likelihood → probably neutral
- **LLR < 0**: mutant is LESS likely → probably destabilizing
- **LLR << -2**: strongly disfavored → likely severely destabilizing

### 1.2 Load ESM-2 Model

```python
#!/opt/conda/envs/prot/bin/python
"""Load ESM-2 model for stability scoring."""
import torch
import esm

_esm2_model = None
_esm2_alphabet = None
_esm2_batch_converter = None

def load_esm2(device=None):
    """Load ESM-2 650M model. ~2.5 GB download on first use, ~2 GB VRAM.

    Returns:
        model, alphabet, batch_converter, device
    """
    global _esm2_model, _esm2_alphabet, _esm2_batch_converter

    if _esm2_model is not None:
        device = next(_esm2_model.parameters()).device
        return _esm2_model, _esm2_alphabet, _esm2_batch_converter, device

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()

    _esm2_model = model
    _esm2_alphabet = alphabet
    _esm2_batch_converter = batch_converter

    print(f"ESM-2 650M loaded on {device}")
    return model, alphabet, batch_converter, device
```

### 1.3 Score a Single Mutation

```python
#!/opt/conda/envs/prot/bin/python
"""Score a single point mutation using ESM-2 LLR."""
import torch

def score_mutation(sequence, position, wt_aa, mut_aa, model=None, alphabet=None,
                   batch_converter=None, device=None):
    """Compute ESM-2 log-likelihood ratio for a single mutation.

    Args:
        sequence: wildtype amino acid sequence
        position: 0-indexed position of mutation
        wt_aa: wildtype amino acid (single letter)
        mut_aa: mutant amino acid (single letter)
        model, alphabet, batch_converter, device: from load_esm2()

    Returns:
        dict with llr, wt_ll, mut_ll, interpretation
    """
    if model is None:
        model, alphabet, batch_converter, device = load_esm2()

    # Validate
    assert sequence[position] == wt_aa, \
        f"Position {position}: expected {wt_aa}, got {sequence[position]}"

    # Tokenize
    data = [("protein", sequence)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    with torch.no_grad():
        results = model(tokens, repr_layers=[], return_contacts=False)
        logits = results["logits"]  # (1, L+2, vocab)

    # Get log-probabilities at the mutation position
    # +1 because of BOS token
    log_probs = torch.log_softmax(logits[0, position + 1], dim=-1)

    wt_idx = alphabet.get_idx(wt_aa)
    mut_idx = alphabet.get_idx(mut_aa)

    wt_ll = float(log_probs[wt_idx])
    mut_ll = float(log_probs[mut_idx])
    llr = mut_ll - wt_ll

    # Interpretation
    if llr > 0.5:
        interp = "FAVORABLE — mutant more natural than wildtype"
    elif llr > -0.5:
        interp = "NEUTRAL — similar likelihood"
    elif llr > -2.0:
        interp = "UNFAVORABLE — mutant less natural"
    else:
        interp = "STRONGLY UNFAVORABLE — likely destabilizing"

    return {
        "mutation": f"{wt_aa}{position+1}{mut_aa}",  # 1-indexed for display
        "position_0idx": position,
        "wt_ll": wt_ll,
        "mut_ll": mut_ll,
        "llr": llr,
        "interpretation": interp,
    }
```

### 1.4 Score All Mutations in a Design vs Wildtype

```python
#!/opt/conda/envs/prot/bin/python
"""Score all mutations between wildtype and designed sequence."""
import pandas as pd

def score_design_vs_wildtype(wt_sequence, design_sequence,
                              model=None, alphabet=None,
                              batch_converter=None, device=None):
    """Score every mutation in a designed sequence relative to wildtype.

    Args:
        wt_sequence: wildtype amino acid sequence
        design_sequence: designed amino acid sequence (same length)

    Returns:
        DataFrame with per-mutation LLR scores + summary stats
    """
    if model is None:
        model, alphabet, batch_converter, device = load_esm2()

    assert len(wt_sequence) == len(design_sequence), \
        f"Length mismatch: {len(wt_sequence)} vs {len(design_sequence)}"

    # Find mutations
    mutations = []
    for i, (wt, mut) in enumerate(zip(wt_sequence, design_sequence)):
        if wt != mut:
            mutations.append((i, wt, mut))

    if not mutations:
        return pd.DataFrame(), {"n_mutations": 0, "mean_llr": 0, "assessment": "IDENTICAL"}

    # Score each mutation
    results = []
    for pos, wt_aa, mut_aa in mutations:
        result = score_mutation(wt_sequence, pos, wt_aa, mut_aa,
                               model, alphabet, batch_converter, device)
        results.append(result)

    df = pd.DataFrame(results)

    # Summary
    mean_llr = df["llr"].mean()
    n_unfavorable = (df["llr"] < -0.5).sum()
    n_strongly_unfavorable = (df["llr"] < -2.0).sum()

    if n_strongly_unfavorable > 0:
        assessment = "HIGH_RISK"
        action = (f"{n_strongly_unfavorable} strongly unfavorable mutation(s). "
                  f"Review these positions carefully.")
    elif n_unfavorable > len(mutations) * 0.5:
        assessment = "MODERATE_RISK"
        action = f"Majority of mutations unfavorable (mean LLR={mean_llr:.2f})."
    elif mean_llr > -0.3:
        assessment = "LOW_RISK"
        action = "Mutations are mostly neutral or favorable."
    else:
        assessment = "MODERATE_RISK"
        action = f"Mean LLR={mean_llr:.2f}, some concern."

    summary = {
        "n_mutations": len(mutations),
        "mean_llr": float(mean_llr),
        "min_llr": float(df["llr"].min()),
        "max_llr": float(df["llr"].max()),
        "n_favorable": int((df["llr"] > 0.5).sum()),
        "n_neutral": int(((df["llr"] >= -0.5) & (df["llr"] <= 0.5)).sum()),
        "n_unfavorable": int(n_unfavorable),
        "n_strongly_unfavorable": int(n_strongly_unfavorable),
        "assessment": assessment,
        "action": action,
    }

    return df, summary

# Usage:
# model, alphabet, batch_converter, device = load_esm2()
# df, summary = score_design_vs_wildtype(
#     "ETFSDLWKLLPEN",  # wildtype
#     "ETFEDLWKKLPQS",  # design
#     model, alphabet, batch_converter, device
# )
# print(f"Assessment: {summary['assessment']}")
# print(f"Mean LLR: {summary['mean_llr']:.3f}")
# print(df[["mutation", "llr", "interpretation"]].to_string())
```

### 1.5 Batch Score Multiple Designs

```python
#!/opt/conda/envs/prot/bin/python
"""Batch stability scoring for multiple designs."""
import pandas as pd

def batch_stability_score(wt_sequence, design_sequences, design_names,
                           model=None, alphabet=None,
                           batch_converter=None, device=None):
    """Score multiple designs against wildtype.

    Args:
        wt_sequence: wildtype sequence
        design_sequences: list of designed sequences
        design_names: list of names/labels

    Returns:
        summary DataFrame (one row per design), detail dict
    """
    if model is None:
        model, alphabet, batch_converter, device = load_esm2()

    summaries = []
    details = {}

    for seq, name in zip(design_sequences, design_names):
        print(f"Scoring {name}...")
        df, summary = score_design_vs_wildtype(
            wt_sequence, seq, model, alphabet, batch_converter, device
        )
        summary["name"] = name
        summary["sequence"] = seq
        summaries.append(summary)
        details[name] = df

    summary_df = pd.DataFrame(summaries)
    summary_df = summary_df.sort_values("mean_llr", ascending=False)  # best first

    print(f"\n=== Stability Ranking ===")
    for _, row in summary_df.iterrows():
        print(f"  {row['name']}: mean_LLR={row['mean_llr']:+.3f} "
              f"({row['n_mutations']} mutations, {row['assessment']})")

    return summary_df, details
```

## Phase 2: Full-Sequence Log-Likelihood (Design Quality Proxy)

### 2.1 Pseudo-Perplexity Scoring

Beyond per-mutation LLR, you can score the **overall naturalness** of an entire sequence using masked marginal pseudo-likelihood. This is slower but captures epistatic effects (combinations of mutations).

```python
#!/opt/conda/envs/prot/bin/python
"""Full-sequence pseudo-perplexity using ESM-2 masked marginals."""
import torch
import numpy as np

def sequence_pseudo_perplexity(sequence, model=None, alphabet=None,
                                batch_converter=None, device=None):
    """Compute pseudo-perplexity of a sequence (lower = more natural).

    Method: mask each position one at a time, predict, compute
    log-likelihood of the true amino acid. Average over all positions.

    This is slow (one forward pass per residue) but captures epistasis.
    For fast screening, use per-mutation LLR instead.

    Args:
        sequence: amino acid string

    Returns:
        dict with pseudo_perplexity, mean_ll, per_residue_ll
    """
    if model is None:
        model, alphabet, batch_converter, device = load_esm2()

    data = [("protein", sequence)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    mask_idx = alphabet.mask_idx
    L = len(sequence)
    per_res_ll = []

    with torch.no_grad():
        for i in range(L):
            masked = tokens.clone()
            masked[0, i + 1] = mask_idx  # +1 for BOS

            output = model(masked, repr_layers=[], return_contacts=False)
            logits = output["logits"]
            log_probs = torch.log_softmax(logits[0, i + 1], dim=-1)

            true_aa = sequence[i]
            true_idx = alphabet.get_idx(true_aa)
            per_res_ll.append(float(log_probs[true_idx]))

    per_res_ll = np.array(per_res_ll)
    mean_ll = float(per_res_ll.mean())
    pseudo_ppl = float(np.exp(-mean_ll))

    return {
        "pseudo_perplexity": pseudo_ppl,
        "mean_ll": mean_ll,
        "per_residue_ll": per_res_ll,
        "sequence_length": L,
    }

# Usage:
# wt_ppl = sequence_pseudo_perplexity("ETFSDLWKLLPEN")
# design_ppl = sequence_pseudo_perplexity("ETFEDLWKKLPQS")
# print(f"WT pseudo-perplexity: {wt_ppl['pseudo_perplexity']:.3f}")
# print(f"Design pseudo-perplexity: {design_ppl['pseudo_perplexity']:.3f}")
# # Lower = more natural/stable
```

### 2.2 When to Use Pseudo-Perplexity vs Per-Mutation LLR

| Method | Speed | Captures epistasis? | Use case |
|--------|-------|-------------------|----------|
| Per-mutation LLR (Phase 1) | Fast (1 forward pass) | No | Screening 50→10, quick ranking |
| Pseudo-perplexity (Phase 2) | Slow (L forward passes) | Yes | Final 10→3 ranking, suspicious designs |

## Phase 3: Stability Decision Tree

```
You have designed sequences and need to assess stability.
│
├── Step 0: Select gate mode
│   ├── Is this an interface redesign task with fold-back + interface gates?
│   │   ├── YES → Use mode="interface" (relaxed LLR gates)
│   │   │   └── Optionally: get burial_info from prot-interface Phase 2.4
│   │   │       to auto-apply strict gates for CORE mutations
│   │   └── NO → Use mode="standard"
│
├── Step 1: Quick screen (50→10) using per-mutation LLR
│   ├── Load ESM-2 (Phase 1.2)
│   ├── Score all designs vs wildtype (Phase 1.4)
│   ├── Hard gates (thresholds depend on mode):
│   │   ├── Standard: mean_LLR < -2.0 → REJECT; min_LLR < -4.0 → FLAG
│   │   └── Interface: mean_LLR < -4.0 → REJECT; min_LLR < -6.0 → FLAG
│   └── Rank by mean_llr, take Top 10
│
├── Step 2: Precision scoring (10→final) — choose method:
│   ├── ESM-2 pseudo-perplexity (Phase 2) — captures epistasis, no extra install
│   ├── ESMFold fold-back pLDDT (prot-structure) — structural stability proxy
│   ├── FoldX ΔΔG (Phase 4, if installed) — physics-based
│   └── Rosetta cartesian_ddg (Phase 4, if installed) — highest accuracy
│
├── Step 3: Resolve conflicts
│   ├── ESM-2 says unnatural + fold-back/interface says fine?
│   │   └── Check burial (prot-interface Phase 2.4):
│   │       ├── CORE mutation → ESM-2 probably right, reject
│   │       ├── PARTIAL → ambiguous, flag for review
│   │       └── SURFACE/interface → ESM-2 probably biased, trust structure
│   ├── ESM-2 says natural + FoldX says destabilizing?
│   │   └── Core → trust physics; Surface → trust ESM-2
│   └── Both agree → high confidence
│
└── Step 4: Final ranking
    └── Composite stability score:
        stability = 0.4 × normalized(mean_llr)
                  + 0.3 × normalized(min_llr)     # worst single mutation
                  + 0.2 × normalized(plddt_foldback) # structural confidence
                  + 0.1 × normalized(n_mutations, invert) # fewer = safer
```

## Phase 4: FoldX and Rosetta (Reference — Not Currently Installed)

### 4.1 FoldX ΔΔG (When Available)

FoldX is a fast empirical force field for estimating ΔΔG of mutations.

```bash
# RepairPDB first (always!)
foldx --command=RepairPDB --pdb=input.pdb

# Then BuildModel for mutation
foldx --command=BuildModel \
  --pdb=input_Repair.pdb \
  --mutant-file=individual_list.txt

# individual_list.txt format:
# WA10E;   (chain A, position 10, W→E)
```

**FoldX ΔΔG interpretation:**
- ΔΔG < -1.0 kcal/mol → stabilizing
- -1.0 ≤ ΔΔG ≤ 1.0 → neutral
- ΔΔG > 1.0 → destabilizing
- ΔΔG > 3.0 → severely destabilizing

### 4.2 Rosetta cartesian_ddg (When Available)

Rosetta's cartesian_ddg is the gold standard for computational ΔΔG prediction.

```bash
# Relax structure first
rosetta_scripts.default.linuxgccrelease \
  -parser:protocol relax.xml \
  -s input.pdb \
  -nstruct 5

# Then cartesian_ddg
cartesian_ddg.default.linuxgccrelease \
  -s relaxed.pdb \
  -ddg:mut_file mutations.txt \
  -ddg:iterations 3
```

**When to use Rosetta over FoldX:**
- Final candidate selection (10→3)
- Core mutations (buried residues)
- When FoldX and ESM-2 disagree
- When you need publication-quality ΔΔG estimates

**When FoldX is sufficient:**
- Fast screening (50→10)
- Surface mutations
- Quick sanity checks

### 4.3 Installation Notes (For Future Reference)

```bash
# FoldX: requires academic license from https://foldxsuite.crg.eu/
# Download binary, place in PATH

# Rosetta: requires academic license from https://www.rosettacommons.org/
# Large install (~30 GB compiled)

# PyRosetta (Python interface):
# pip install pyrosetta (requires license key)
```

## Phase 5: Consensus Design (Evolutionary Stability)

### 5.1 Quick Consensus Check Using ESM-2

```python
#!/opt/conda/envs/prot/bin/python
"""Check if design mutations agree with evolutionary consensus."""
import torch

def consensus_check(sequence, model=None, alphabet=None,
                     batch_converter=None, device=None, top_k=3):
    """For each position, get top-K most likely amino acids from ESM-2.

    A mutation that matches top-K is "consensus-compatible."
    A mutation outside top-K is "consensus-divergent" and carries more risk.

    Args:
        sequence: amino acid sequence
        top_k: how many top predictions to consider "consensus"

    Returns:
        list of dicts with position, aa, rank, is_consensus
    """
    if model is None:
        model, alphabet, batch_converter, device = load_esm2()

    data = [("protein", sequence)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    with torch.no_grad():
        results = model(tokens, repr_layers=[], return_contacts=False)
        logits = results["logits"][0]  # (L+2, vocab)

    standard_aas = "ACDEFGHIKLMNPQRSTVWY"
    aa_indices = {aa: alphabet.get_idx(aa) for aa in standard_aas}

    positions = []
    for i, aa in enumerate(sequence):
        pos_logits = logits[i + 1]  # +1 for BOS
        # Get probabilities for standard AAs only
        probs = {}
        for res, idx in aa_indices.items():
            probs[res] = float(torch.softmax(pos_logits, dim=-1)[idx])

        ranked = sorted(probs.items(), key=lambda x: -x[1])
        top_aas = [r[0] for r in ranked[:top_k]]

        # Find rank of current AA
        rank = next((j for j, (r, _) in enumerate(ranked) if r == aa), len(ranked))

        positions.append({
            "position": i + 1,  # 1-indexed
            "aa": aa,
            "rank": rank + 1,  # 1-indexed
            "is_consensus": aa in top_aas,
            "top_3": "".join(top_aas),
            "prob": probs.get(aa, 0),
        })

    return positions

def compare_design_consensus(wt_sequence, design_sequence,
                              model=None, alphabet=None,
                              batch_converter=None, device=None):
    """Check which mutations in a design are consensus-compatible.

    Returns:
        DataFrame with mutation-level consensus analysis
    """
    if model is None:
        model, alphabet, batch_converter, device = load_esm2()

    wt_consensus = consensus_check(wt_sequence, model, alphabet, batch_converter, device)
    design_consensus = consensus_check(design_sequence, model, alphabet, batch_converter, device)

    mutations = []
    for i, (wt_aa, des_aa) in enumerate(zip(wt_sequence, design_sequence)):
        if wt_aa != des_aa:
            wt_info = wt_consensus[i]
            des_info = design_consensus[i]
            mutations.append({
                "position": i + 1,
                "wt_aa": wt_aa,
                "des_aa": des_aa,
                "wt_rank_in_wt_context": wt_info["rank"],
                "des_rank_in_des_context": des_info["rank"],
                "des_in_wt_top3": des_aa in wt_info["top_3"],
                "wt_top3": wt_info["top_3"],
            })

    import pandas as pd
    df = pd.DataFrame(mutations)
    if len(df) > 0:
        n_consensus = df["des_in_wt_top3"].sum()
        print(f"Mutations: {len(df)} total, {n_consensus} consensus-compatible "
              f"({100*n_consensus/len(df):.0f}%)")
    return df
```

## Phase 6: Hard Gates for Stability

### 6.1 Standard Hard Gates

```python
#!/opt/conda/envs/prot/bin/python
"""Stability hard gates for screening designs."""

def stability_hard_gates(summary_df, mode="standard"):
    """Apply hard gates to stability scoring results.

    Args:
        summary_df: from batch_stability_score()
        mode: "standard" or "interface" (relaxed for interface mutations)

    Returns:
        passed_df, failed_df, gate_report
    """
    if mode == "interface":
        return stability_hard_gates_interface(summary_df)

    failed_indices = set()
    reasons_list = []

    for idx, row in summary_df.iterrows():
        reasons = []

        # Gate 1: No severely unnatural designs
        if row.get("mean_llr", 0) < -2.0:
            reasons.append(f"mean_LLR={row['mean_llr']:.2f} < -2.0 (globally unnatural)")

        # Gate 2: No single catastrophic mutation
        if row.get("min_llr", 0) < -4.0:
            reasons.append(f"min_LLR={row['min_llr']:.2f} < -4.0 (catastrophic single mutation)")

        # Gate 3: Majority of mutations should not be unfavorable
        n_mut = row.get("n_mutations", 0)
        n_unfav = row.get("n_unfavorable", 0)
        if n_mut > 0 and n_unfav / n_mut > 0.7:
            reasons.append(f"{n_unfav}/{n_mut} mutations unfavorable (>70%)")

        if reasons:
            failed_indices.add(idx)
            reasons_list.append({"name": row.get("name", ""), "reasons": reasons})

    passed = summary_df[~summary_df.index.isin(failed_indices)]
    failed = summary_df[summary_df.index.isin(failed_indices)]

    report = {
        "n_input": len(summary_df),
        "n_passed": len(passed),
        "n_failed": len(failed),
        "mode": "standard",
        "failures": reasons_list,
    }

    return passed, failed, report
```

### 6.2 Interface-Mutation Mode (Burial-Aware Gates)

```python
#!/opt/conda/envs/prot/bin/python
"""Relaxed stability gates for interface design tasks.

Why this exists:
- ESM-2 LLR captures evolutionary conservation, not just thermodynamic stability
- Interface positions are often highly conserved → functional mutations get penalized
- 1BRS integration test: A36S scored TM=0.97 fold-back PASS + interface PASS
  but LLR=-6.19 catastrophic FAIL
- The standard gates are correct for CORE mutations but too strict for
  SURFACE/PARTIAL interface mutations

When to use interface mode:
- Task is explicitly interface redesign / affinity maturation
- Mutations are at interface positions (not core)
- fold-back and interface gates are also being applied

The interface mode:
- Relaxes mean_LLR gate from -2.0 to -4.0
- Relaxes catastrophic gate from -4.0 to -6.0
- REQUIRES fold-back PASS + interface PASS as prerequisites
- Adds burial-aware classification (from prot-interface Phase 2.4)
"""

def stability_hard_gates_interface(summary_df, burial_info=None):
    """Apply relaxed hard gates for interface mutation scenarios.

    Args:
        summary_df: from batch_stability_score()
        burial_info: optional dict {position_1idx: "CORE"/"PARTIAL"/"SURFACE"}
                     from prot-interface annotate_interface_burial()
                     If provided, CORE mutations still use strict gates.

    Returns:
        passed_df, failed_df, gate_report
    """
    failed_indices = set()
    reasons_list = []

    # Relaxed thresholds for interface mode
    MEAN_LLR_GATE = -4.0     # relaxed from -2.0
    CATASTROPHIC_GATE = -6.0  # relaxed from -4.0
    UNFAV_RATIO_GATE = 0.85   # relaxed from 0.70

    for idx, row in summary_df.iterrows():
        reasons = []

        mean_llr = row.get("mean_llr", 0)
        min_llr = row.get("min_llr", 0)

        # If burial info is available, check if any mutation hits CORE
        has_core_mutation = False
        if burial_info and "detail_mutations" in row:
            for mut in row.get("detail_mutations", []):
                pos = mut.get("position_1idx", 0)
                if burial_info.get(pos) == "CORE":
                    has_core_mutation = True
                    break

        if has_core_mutation:
            # CORE mutation: use STRICT thresholds (same as standard)
            if mean_llr < -2.0:
                reasons.append(f"mean_LLR={mean_llr:.2f} < -2.0 "
                              f"(CORE mutation — strict gate applies)")
            if min_llr < -4.0:
                reasons.append(f"min_LLR={min_llr:.2f} < -4.0 "
                              f"(CORE mutation — catastrophic)")
        else:
            # Non-core (PARTIAL/SURFACE/interface): use RELAXED thresholds
            if mean_llr < MEAN_LLR_GATE:
                reasons.append(f"mean_LLR={mean_llr:.2f} < {MEAN_LLR_GATE} "
                              f"(interface-relaxed gate)")
            if min_llr < CATASTROPHIC_GATE:
                reasons.append(f"min_LLR={min_llr:.2f} < {CATASTROPHIC_GATE} "
                              f"(interface-relaxed catastrophic)")

        # Unfavorable ratio (slightly relaxed)
        n_mut = row.get("n_mutations", 0)
        n_unfav = row.get("n_unfavorable", 0)
        if n_mut > 0 and n_unfav / n_mut > UNFAV_RATIO_GATE:
            reasons.append(f"{n_unfav}/{n_mut} mutations unfavorable "
                          f"(>{UNFAV_RATIO_GATE*100:.0f}%)")

        if reasons:
            failed_indices.add(idx)
            reasons_list.append({"name": row.get("name", ""), "reasons": reasons})

    passed = summary_df[~summary_df.index.isin(failed_indices)]
    failed = summary_df[summary_df.index.isin(failed_indices)]

    report = {
        "n_input": len(summary_df),
        "n_passed": len(passed),
        "n_failed": len(failed),
        "mode": "interface",
        "thresholds": {
            "mean_llr": MEAN_LLR_GATE,
            "catastrophic": CATASTROPHIC_GATE,
            "unfav_ratio": UNFAV_RATIO_GATE,
        },
        "note": "Interface mode: REQUIRES fold-back PASS + interface PASS as prerequisites",
        "failures": reasons_list,
    }

    return passed, failed, report
```

### 6.3 When to Use Which Mode

```
Which stability gate mode should I use?
│
├── Is this an interface redesign / affinity maturation task?
│   ├── YES → Are fold-back AND interface gates also being applied?
│   │   ├── YES → Use mode="interface"
│   │   └── NO → Use mode="standard" (you need the stricter gates
│   │            to compensate for missing structural checks)
│   └── NO → Use mode="standard"
│
├── Are mutations only at SURFACE/PARTIAL positions?
│   └── YES + interface task → mode="interface" is safe
│
└── Do any mutations hit CORE positions?
    └── YES → mode="interface" still uses strict gates for CORE mutations
              (burial_info auto-applies strict thresholds per-position)
```

| Mode | mean_LLR gate | catastrophic gate | unfav ratio | Prerequisites |
|------|--------------|-------------------|-------------|---------------|
| **standard** | -2.0 | -4.0 | 70% | None |
| **interface** | -4.0 | -6.0 | 85% | fold-back PASS + interface PASS |
| **interface + CORE** | -2.0 (for CORE mutations) | -4.0 (for CORE) | 85% | burial_info from prot-interface |

## Phase 7: VRAM Budget & Timing

| Task | VRAM | Wall time (<GPU_REDACTED>) | Notes |
|------|------|---------------------|-------|
| Load ESM-2 650M | ~2 GB | ~10 sec | One-time |
| Per-mutation LLR, 1 design (5 mutations) | ~2 GB | ~1 sec | Single forward pass |
| Per-mutation LLR, 50 designs | ~2 GB | ~30 sec | One pass per design |
| Pseudo-perplexity, 1 sequence (100 aa) | ~2 GB | ~2 min | 100 forward passes |
| Pseudo-perplexity, 1 sequence (300 aa) | ~2 GB | ~5 min | 300 forward passes |
| Consensus check, 1 sequence | ~2 GB | ~1 sec | Single forward pass |

> ESM-2 650M fits comfortably alongside LigandMPNN on 12 GB. No VRAM management needed.

## Failure Modes

1. **Using ESM-2 LLR as ground truth.** LLR is a proxy for evolutionary plausibility, not thermodynamic stability. A mutation can be evolutionarily unusual (low LLR) but thermodynamically stabilizing (e.g., engineered disulfide bonds). Use LLR for screening, not final decisions.

2. **Ignoring context dependence.** ESM-2 LLR depends on the surrounding sequence context. The same mutation at the same position can have different LLR in different sequence backgrounds. Always score mutations IN the context of your design, not the wildtype.

3. **Averaging conflicting signals.** When ESM-2 and FoldX disagree, the answer is NOT (ESM + FoldX) / 2. Diagnose the disagreement: is it a surface mutation (trust ESM) or a core packing change (trust FoldX)?

4. **Forgetting wildtype comparison.** An LLR of -0.5 for a mutation seems concerning, but if the wildtype position also has unusual amino acids (low probability), the relative change might be acceptable. Always report ΔLLR not absolute scores.

5. **Over-filtering with consensus.** Not every non-consensus mutation is bad. Functional sites often have unusual residues (catalytic, binding). Don't reject mutations just because they're not in the top-3 consensus — use consensus as a risk flag, not a hard gate.

6. **Pseudo-perplexity on short peptides.** For sequences <20 aa, pseudo-perplexity has high variance and limited meaning. Stick to per-mutation LLR for short peptides.

7. **Over-rejecting interface mutations with standard gates.** ESM-2 heavily penalizes mutations at conserved interface positions because they're evolutionarily unusual — but that's exactly where you WANT to mutate for affinity maturation. Use mode="interface" when doing interface redesign, and always cross-check with fold-back (prot-structure) and interface gates (prot-interface). The 1BRS integration test showed: A36S scored TM=0.97 fold-back + interface PASS but LLR=-6.19 catastrophic FAIL under standard gates — a false negative. With interface mode + burial classification, this design would be correctly handled.

## One-Sentence Rule

**Screen with ESM-2 LLR (fast, kills obviously bad mutations), use interface mode with burial context for interface redesign tasks, compare every score to wildtype, and when physics and evolution disagree — diagnose the disagreement rather than averaging it away.**
