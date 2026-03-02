---
name: prot-seqdesign
description: Protein sequence design using LigandMPNN and ESM-IF. Provides complete CLI reference, input preparation, output parsing, temperature/diversity tuning, multi-chain and interface-only design, and quality gates. Triggers whenever you need to design sequences for a given protein backbone — fixed-backbone design, interface redesign, or stability optimization via sequence.
homepage: https://github.com/dauparas/LigandMPNN
metadata: { "openclaw": { "emoji": "🧬", "requires": { "bins": ["python3"], "python": ["biopython", "fair-esm", "torch", "numpy", "pandas", "prody", "torch_geometric", "tabulate"] } } }
---

# Protein Sequence Design: LigandMPNN + ESM-IF

Sequence design is the core generative step in protein engineering. Given a backbone structure, produce amino acid sequences that fold into that structure, maintain stability, and satisfy functional constraints (binding, catalysis, solubility). This skill makes you competent at running the tools — you already understand the methodology.

## When to Use

- Designing sequences for a **fixed backbone** (from PDB, ESMFold, or RFdiffusion)
- **Interface redesign** — only mutating residues at a binding interface
- **Stability optimization** — redesigning core or surface residues
- **Consensus-guided design** — combining computational design with evolutionary information
- Generating **diverse sequence libraries** for experimental screening
- **Cross-validating** designs between LigandMPNN and ESM-IF

## Core Philosophy

1. **Structure quality gates before sequence design.** Never design on a bad backbone. Check pLDDT/PAE/clashes BEFORE running LigandMPNN. Garbage in, garbage out — but with confident-looking sequences.
2. **Compare to wildtype, always.** Every design metric (ESM log-likelihood, predicted stability, folding confidence) must be compared to the wildtype sequence on the same backbone. A design that scores -3.2 means nothing without knowing wildtype scores -3.5.
3. **Temperature is your diversity knob, not your quality knob.** Low temperature (0.1) gives high recovery but low diversity. High temperature (0.3+) gives diversity but lower per-sequence confidence. The right temperature depends on your goal, not on "what's best."
4. **Design the minimum region necessary.** Don't redesign the entire protein when you only need to optimize the interface. Every unnecessary mutation is a risk with no upside.
5. **Two tools are better than one.** LigandMPNN and ESM-IF use different architectures. Agreement between them is a confidence signal; disagreement is an uncertainty flag.

## Phase 0: Environment Verification

**Run this BEFORE any design task.** This phase detects which LigandMPNN version/fork you have and verifies all dependencies.

```python
#!/opt/conda/envs/prot/bin/python
"""Phase 0: Environment verification for prot-seqdesign skill.

Run this once per session or after any package update.
It detects your LigandMPNN version and sets the correct CLI flags.
"""
import subprocess, sys, os, importlib

PYTHON = "/opt/conda/envs/prot/bin/python"
MPNN_DIR = "<ABS_PATH_REDACTED>"
# If you vendor'd a local copy, override:
# MPNN_DIR = "<ABS_PATH_REDACTED>"

print("=== Phase 0: Environment Verification ===\n")

# 1. Core dependencies
deps = {
    "torch": None, "numpy": None, "pandas": None,
    "Bio.PDB": "biopython", "esm": "fair-esm",
    "prody": None, "tabulate": None,
}
missing = []
for mod, pkg in deps.items():
    try:
        importlib.import_module(mod.split(".")[0])
    except ImportError:
        missing.append(pkg or mod)
if missing:
    print(f"❌ Missing packages: {missing}")
    print(f"   Fix: /opt/conda/envs/prot/bin/pip install {' '.join(missing)}")
else:
    print("✅ Core dependencies OK")

# 2. PyG dependencies (for ESM-IF)
pyg_deps = ["torch_geometric", "torch_scatter", "torch_sparse", "torch_cluster"]
pyg_missing = []
for mod in pyg_deps:
    try:
        importlib.import_module(mod)
    except ImportError:
        pyg_missing.append(mod)
if pyg_missing:
    print(f"⚠️  PyG dependencies missing: {pyg_missing}")
    print(f"   Fix: /opt/conda/envs/prot/bin/pip install torch_geometric")
    print(f"   Then: /opt/conda/envs/prot/bin/pip install torch_scatter torch_sparse "
          f"torch_cluster torch_spline_conv "
          f"-f https://data.pyg.org/whl/torch-2.6.0+cu124.html")
else:
    print("✅ PyG dependencies OK (ESM-IF ready)")

# 3. GPU
import torch
if torch.cuda.is_available():
    print(f"✅ GPU: {torch.cuda.get_device_name(0)} "
          f"({torch.cuda.get_device_properties(0).total_memory/1024**3:.1f} GB)")
else:
    print("⚠️  No GPU — will be slow")

# 4. LigandMPNN version detection
run_py = os.path.join(MPNN_DIR, "run.py")
if not os.path.isfile(run_py):
    print(f"❌ LigandMPNN not found at {run_py}")
else:
    # Detect CLI flag style by parsing --help
    result = subprocess.run(
        [PYTHON, run_py, "-h"],
        capture_output=True, text=True, timeout=30
    )
    help_text = result.stdout + result.stderr

    if "--number_of_batches" in help_text:
        cli_style = "NEW"
        print("✅ LigandMPNN detected: NEW-style CLI")
        print("   Use: --number_of_batches N --batch_size 1 --temperature T")
    elif "--num_seq_per_target" in help_text:
        cli_style = "LEGACY"
        print("✅ LigandMPNN detected: LEGACY-style CLI")
        print("   Use: --num_seq_per_target N --sampling_temp T")
    else:
        cli_style = "UNKNOWN"
        print("⚠️  LigandMPNN CLI style unknown. Run: python run.py -h")

    # Check weights
    weights_dir = os.path.join(MPNN_DIR, "model_params")
    if os.path.isdir(weights_dir):
        n = len(os.listdir(weights_dir))
        print(f"✅ Model weights: {n} files in {weights_dir}")
    else:
        print(f"❌ No model_params/ directory in {MPNN_DIR}")

# 5. biotite compatibility check (for ESM-IF)
try:
    from biotite.structure import filter_backbone
    print("✅ biotite filter_backbone OK")
except ImportError:
    try:
        from biotite.structure import filter_peptide_backbone
        print("⚠️  biotite uses filter_peptide_backbone (needs runtime patch for ESM-IF)")
        print("   The ESM-IF scoring function in this skill handles this automatically.")
    except ImportError:
        print("❌ biotite backbone filter not found")

print("\n=== Phase 0 complete ===")
```

### Phase 0 Output Variables

After running Phase 0, set these for the rest of the pipeline:

```python
# === CONFIGURE BASED ON PHASE 0 OUTPUT ===
PYTHON = "/opt/conda/envs/prot/bin/python"
MPNN_DIR = "<ABS_PATH_REDACTED>"
MPNN_RUN = f"{MPNN_DIR}/run.py"
MPNN_WEIGHTS = f"{MPNN_DIR}/model_params/proteinmpnn_v_48_020.pt"

# Set by Phase 0 detection:
CLI_STYLE = "NEW"  # or "LEGACY" — determines which flags to use
```

## Installation & Paths

### Paths in this container

```
Python:         /opt/conda/envs/prot/bin/python
LigandMPNN:     <ABS_PATH_REDACTED>
MPNN weights:   <ABS_PATH_REDACTED>
ESM-IF model:   auto-downloaded on first import (~500 MB)
```

### Dependency installation (if Phase 0 reports missing)

```bash
# Core
/opt/conda/envs/prot/bin/pip install biopython biotite fair-esm prody tabulate tmtools

# PyTorch Geometric (required for ESM-IF)
/opt/conda/envs/prot/bin/pip install torch_geometric
/opt/conda/envs/prot/bin/pip install \
  torch_scatter torch_sparse torch_cluster torch_spline_conv \
  -f https://data.pyg.org/whl/torch-2.6.0+cu124.html

# If biotite version mismatch — the ESM-IF function below handles this at runtime
```

## Phase 1: Input Preparation (PDB → Clean Structure)

### 1.1 Download and Clean a PDB

```python
#!/opt/conda/envs/prot/bin/python
"""Download PDB, remove water/heteroatoms, select chains."""
import os, requests
from Bio.PDB import PDBParser, PDBIO, Select

def download_pdb(pdb_id, out_path):
    """Download PDB from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with open(out_path, "w") as f:
        f.write(r.text)
    print(f"Downloaded {pdb_id} → {out_path}")

class CleanSelect(Select):
    """Keep only specified chains, remove water and heteroatoms."""
    def __init__(self, keep_chains=None):
        self.keep_chains = set(keep_chains) if keep_chains else None
    def accept_chain(self, chain):
        if self.keep_chains is None:
            return True
        return chain.id in self.keep_chains
    def accept_residue(self, residue):
        hetflag = residue.id[0]
        if residue.get_resname() == "HOH":
            return False
        if hetflag != " ":
            return False
        return True

def clean_pdb(in_pdb, out_pdb, keep_chains=None):
    """Clean PDB: remove water, heteroatoms, select chains."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", in_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb, CleanSelect(keep_chains))
    print(f"Cleaned → {out_pdb} (chains: {keep_chains or 'all'})")
```

### 1.2 Extract Interface Residues

```python
#!/opt/conda/envs/prot/bin/python
"""Identify interface residues between two chains."""
import numpy as np
from Bio.PDB import PDBParser

def get_interface_residues(pdb_path, chain_a, chain_b, cutoff=5.0):
    """Find residues in chain_a within cutoff Å of chain_b (heavy atoms).

    Args:
        pdb_path: path to cleaned PDB
        chain_a: chain ID to find interface residues on
        chain_b: partner chain ID
        cutoff: distance threshold in Ångströms (default 5.0)

    Returns:
        list of (chain_id, resseq) tuples for interface residues in chain_a
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    atoms_b = [a for a in model[chain_b].get_atoms() if a.element != "H"]
    coords_b = np.array([a.coord for a in atoms_b])

    interface = []
    for res in model[chain_a].get_residues():
        if res.id[0] != " ":
            continue
        atoms_a = [a for a in res.get_atoms() if a.element != "H"]
        if not atoms_a:
            continue
        coords_a = np.array([a.coord for a in atoms_a])
        dists = np.sqrt(((coords_a[:, None] - coords_b[None, :]) ** 2).sum(-1))
        if dists.min() <= cutoff:
            interface.append((chain_a, res.id[1]))

    print(f"Interface residues on chain {chain_a}: {len(interface)} "
          f"(cutoff={cutoff}Å)")
    return interface

def format_residues_for_mpnn(interface_residues):
    """Convert [(chain, resseq), ...] to LigandMPNN --redesigned_residues format.

    Returns:
        str like "B17 B18 B19 ..."
    """
    return " ".join(f"{chain}{resseq}" for chain, resseq in interface_residues)
```

## Phase 2: LigandMPNN — Complete CLI Reference

### ⚠️ CRITICAL: Version-Dependent CLI Flags

LigandMPNN has evolved. Different versions use different flag names. **Always run Phase 0 first** to detect your version, or run `python run.py -h`.

| Concept | NEW-style flags | LEGACY-style flags |
|---------|-----------------|-------------------|
| Number of sequences | `--number_of_batches N --batch_size 1` | `--num_seq_per_target N` |
| Temperature | `--temperature "0.1"` | `--sampling_temp "0.1"` |
| Output score field | `overall_confidence` (higher=better) | `score` (lower=better) |
| Recovery field | `seq_rec` | `seq_recovery` |
| Redesigned residues | `--redesigned_residues "B17 B18..."` | `--redesigned_residues "B17 B18..."` (same) |
| Fixed residues | `--fixed_residues "A10 A20..."` | `--fixed_residues "A10 A20..."` (same) |
| Tied positions | `--tied_positions "A1 B1, A2 B2..."` | `--tied_positions "A1 B1, A2 B2..."` (same) |

> **Score direction differs!** NEW: `overall_confidence` higher = better. LEGACY: `score` lower = better. Get this wrong and your ranking inverts.

### 2.1 Run Commands (Both Styles)

#### NEW-style CLI

```bash
/opt/conda/envs/prot/bin/python <ABS_PATH_REDACTED> \
  --model_type "protein_mpnn" \
  --checkpoint_protein_mpnn "<ABS_PATH_REDACTED>" \
  --pdb_path "input_clean.pdb" \
  --out_folder "output/" \
  --number_of_batches 50 \
  --batch_size 1 \
  --temperature "0.1" \
  --seed 42
```

#### LEGACY-style CLI

```bash
/opt/conda/envs/prot/bin/python <ABS_PATH_REDACTED> \
  --model_type "protein_mpnn" \
  --checkpoint_protein_mpnn "<ABS_PATH_REDACTED>" \
  --pdb_path "input_clean.pdb" \
  --out_folder "output/" \
  --num_seq_per_target 50 \
  --sampling_temp "0.1" \
  --seed 42 \
  --batch_size 1
```

### 2.2 Key Parameters

| Parameter | What it does | Recommended |
|-----------|-------------|-------------|
| `--model_type` | Model variant | `"protein_mpnn"` standard; `"ligand_mpnn"` if ligand present |
| `--checkpoint_protein_mpnn` | Weight file | `proteinmpnn_v_48_020.pt` (noise=0.20, good general use) |
| `--pdb_path` | Input PDB (single) | Cleaned PDB, no water/HETATM |
| `--out_folder` | Output directory | Will contain seqs/*.fa |
| Sequences | NEW: `--number_of_batches 50 --batch_size 1` / LEGACY: `--num_seq_per_target 50` | 50 explore, 10-20 focused |
| Temperature | NEW: `--temperature "0.1"` / LEGACY: `--sampling_temp "0.1"` | 0.1 conservative, 0.2 balanced, 0.3 diverse |
| `--seed` | Random seed | Set for reproducibility |

### 2.3 Controlling What Gets Designed

#### Design only specific chains

```bash
# Design chain B only; chain A is fixed
... --redesigned_residues "B" ...
```

> `--redesigned_residues "B"` = redesign ALL residues on chain B.

#### Design specific positions only

```bash
# Design only positions 17-29 on chain B
... --redesigned_residues "B17 B18 B19 B20 B22 B23 B25 B26 B27 B28 B29" ...
```

> **Format**: space-separated `<chain><resnum>`. Get residue numbers from Phase 1.2.

#### Fix specific positions (design everything else)

```bash
# Fix catalytic residues; design the rest
... --fixed_residues "A10 A20 A30" ...
```

#### Tied positions (symmetric design)

```bash
# Homo-dimer: tie corresponding positions
... --tied_positions "A1 B1, A2 B2, A3 B3" ...
```

### 2.4 Model Weight Selection Guide

```
Is there a bound ligand/cofactor/metal in the PDB?
├── YES → --model_type "ligand_mpnn"
│         --checkpoint_ligand_mpnn ".../ligandmpnn_v_32_020_25.pt"
└── NO
    ├── Is it a membrane protein?
    │   └── YES → --checkpoint_protein_mpnn ".../global_label_membrane_mpnn_v_48_020.pt"
    └── NO
        ├── Is solubility/expression a priority?
        │   ├── YES → --checkpoint_protein_mpnn ".../solublempnn_v_48_020.pt"
        │   └── NO  → --checkpoint_protein_mpnn ".../proteinmpnn_v_48_020.pt"  (DEFAULT)
        └──
```

### 2.5 Output Format

LigandMPNN writes FASTA files to `out_folder/seqs/`.

**NEW-style output** (FASTA headers):

```
>1YCR_clean, overall_confidence=..., ligand_confidence=..., seq_rec=...
MKTLYILGLVGAIALSMVTVQGKNL...
>T=0.1, sample=1, overall_confidence=0.8765, ligand_confidence=0.0, seq_rec=0.6364
MKTLYILGLVGAIALSMVTTQGKNL...
```

**LEGACY-style output**:

```
>1YCR_clean, score=1.2345, global_score=1.3456, seq_recovery=0.4567
MKTLYILGLVGAIALSMVTVQGKNL...
>T=0.1, sample=1, score=0.9876, global_score=1.0234, seq_recovery=0.4321
MKTLYILGLVGAIALSMVTTQGKNL...
```

- **Line 1 (first record)**: wildtype/input sequence (reference)
- NEW: `overall_confidence` higher = better; LEGACY: `score` lower = better
- `seq_rec`/`seq_recovery`: fraction of designed positions matching wildtype

### 2.6 Parse LigandMPNN Output (Version-Adaptive)

```python
#!/opt/conda/envs/prot/bin/python
"""Parse LigandMPNN FASTA output — handles both NEW and LEGACY formats."""
import re
import pandas as pd

def parse_mpnn_fasta(fasta_path):
    """Parse LigandMPNN output FASTA. Auto-detects NEW vs LEGACY format.

    Returns:
        list of dicts, format_style ("NEW" or "LEGACY")
    """
    records = []
    with open(fasta_path) as f:
        lines = f.readlines()

    header, seq = None, ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append(_parse_record(header, seq))
            header = line[1:]
            seq = ""
        else:
            seq += line
    if header is not None:
        records.append(_parse_record(header, seq))

    # Detect format
    first_header = records[0]["header"] if records else ""
    if "overall_confidence" in first_header:
        style = "NEW"
    elif "score=" in first_header and "overall_confidence" not in first_header:
        style = "LEGACY"
    else:
        style = "UNKNOWN"

    return records, style

def _parse_record(header, sequence):
    """Extract fields from a FASTA header (either format)."""
    rec = {"header": header, "sequence": sequence, "is_wildtype": False}

    # NEW-style fields
    m = re.search(r"overall_confidence=([\d.]+)", header)
    if m:
        rec["overall_confidence"] = float(m.group(1))

    m = re.search(r"ligand_confidence=([\d.]+)", header)
    if m:
        rec["ligand_confidence"] = float(m.group(1))

    # LEGACY-style fields
    m = re.search(r"(?<![a-z_])score=([\d.]+)", header)
    if m:
        rec["score"] = float(m.group(1))

    m = re.search(r"global_score=([\d.]+)", header)
    if m:
        rec["global_score"] = float(m.group(1))

    # Recovery (both styles)
    m = re.search(r"seq_rec(?:overy)?=([\d.]+)", header)
    if m:
        rec["seq_recovery"] = float(m.group(1))

    # Sample index
    m = re.search(r"sample=(\d+)", header)
    if m:
        rec["sample_idx"] = int(m.group(1))
    else:
        rec["is_wildtype"] = True
        rec["sample_idx"] = 0

    return rec

def mpnn_output_to_df(fasta_path):
    """Parse into a sorted DataFrame. Auto-detects format and sorts correctly."""
    records, style = parse_mpnn_fasta(fasta_path)
    df = pd.DataFrame(records)

    if style == "NEW" and "overall_confidence" in df.columns:
        df = df.sort_values("overall_confidence", ascending=False)  # higher = better
        print(f"Format: NEW (overall_confidence, higher=better)")
    elif style == "LEGACY" and "score" in df.columns:
        df = df.sort_values("score", ascending=True)  # lower = better
        print(f"Format: LEGACY (score, lower=better)")
    else:
        print(f"Format: {style} — check sorting manually")

    return df, style

def compute_interface_recovery(wt_seq, design_seq, interface_positions):
    """Compute recovery specifically over interface residues.

    Args:
        wt_seq: wildtype full-chain sequence (str)
        design_seq: designed full-chain sequence (str)
        interface_positions: list of 0-indexed positions in the sequence

    Returns:
        float: fraction of interface positions matching wildtype
    """
    if len(wt_seq) != len(design_seq):
        raise ValueError(f"Sequence length mismatch: {len(wt_seq)} vs {len(design_seq)}")
    matches = sum(1 for i in interface_positions if wt_seq[i] == design_seq[i])
    return matches / len(interface_positions) if interface_positions else 0.0
```

## Phase 3: ESM-IF — Inverse Folding Cross-Validation

### 3.1 Score a Sequence Against a Structure

```python
#!/opt/conda/envs/prot/bin/python
"""Score sequence-structure compatibility using ESM-IF.

Handles biotite API differences and GPU device placement.
Tested on fair-esm 2.0.0 + biotite 1.6.0.
"""
import torch
import esm
import esm.inverse_folding
import warnings

def _patch_biotite_if_needed():
    """Patch biotite filter_backbone → filter_peptide_backbone if needed."""
    import biotite.structure as bs
    if not hasattr(bs, 'filter_backbone') and hasattr(bs, 'filter_peptide_backbone'):
        bs.filter_backbone = bs.filter_peptide_backbone
        warnings.warn("Patched biotite: filter_backbone → filter_peptide_backbone")

def load_esmif_model(device=None):
    """Load ESM-IF model. ~500 MB download on first use.

    Args:
        device: "cuda" or "cpu" (default: auto-detect)

    Returns:
        model, alphabet, device
    """
    _patch_biotite_if_needed()

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().to(device)
    return model, alphabet, device

def score_sequence_on_structure(model, alphabet, pdb_path, chain_id,
                                 sequence=None, device="cuda"):
    """Score a sequence's compatibility with a backbone structure.

    Args:
        model, alphabet: from load_esmif_model()
        pdb_path: path to PDB file
        chain_id: which chain to score
        sequence: amino acid sequence (if None, uses native sequence from PDB)
        device: torch device string

    Returns:
        dict with:
        - log_likelihood: total log-likelihood (higher = better fit)
        - avg_ll: mean log-likelihood per residue
        - sequence: the scored sequence
        - native_sequence: native sequence from PDB
    """
    _patch_biotite_if_needed()

    # Load structure and extract coords
    structure = esm.inverse_folding.util.load_structure(pdb_path, chain_id)
    coords, native_seq = esm.inverse_folding.util.extract_coords_from_structure(structure)

    if sequence is None:
        sequence = native_seq

    # Prepare batch using the batch converter
    batch_converter = esm.inverse_folding.util.CoordBatchConverter(alphabet)
    batch = [(coords, None, sequence)]
    coords_batch, confidence, strs, tokens, padding_mask = batch_converter(batch)

    # Move ALL tensors to the same device as model
    coords_batch = coords_batch.to(device)
    confidence = confidence.to(device)
    tokens = tokens.to(device)
    padding_mask = padding_mask.to(device)

    with torch.no_grad():
        # ESM-IF forward: (coords, padding_mask, confidence, prev_output_tokens)
        logits = model.forward(
            coords_batch,
            padding_mask,
            confidence,
            tokens[:, :-1],  # prev_output_tokens (teacher forcing, shift right)
        )
        # Per-residue log-likelihood
        log_probs = torch.log_softmax(logits, dim=-1)
        target = tokens[:, 1:]  # target tokens (shifted)
        per_res_ll = log_probs.gather(2, target.unsqueeze(-1)).squeeze(-1)
        # Mask padding
        mask = ~padding_mask[:, 1:].to(device)
        per_res_ll = per_res_ll * mask

    per_res = per_res_ll[0].cpu().numpy()
    n_valid = int(mask[0].sum().item())
    total = float(per_res.sum())
    avg = total / n_valid if n_valid > 0 else 0.0

    return {
        "log_likelihood": total,
        "avg_ll": avg,
        "n_residues": n_valid,
        "sequence": sequence,
        "native_sequence": native_seq,
    }

# Usage:
# model, alphabet, device = load_esmif_model()
# wt_result = score_sequence_on_structure(model, alphabet, "complex.pdb", "B", device=device)
# design_result = score_sequence_on_structure(model, alphabet, "complex.pdb", "B",
#                                              sequence="ETFEDLWKKLPQS", device=device)
# delta = design_result["avg_ll"] - wt_result["avg_ll"]
# print(f"WT avg_ll: {wt_result['avg_ll']:.4f}")
# print(f"Design avg_ll: {design_result['avg_ll']:.4f}")
# print(f"Delta: {delta:+.4f} ({'better' if delta > 0 else 'worse'})")
```

### 3.2 ESM-IF VRAM Budget

| Input size | Approx VRAM | 12 GB strategy |
|-----------|-------------|---------------|
| ≤300 aa | ~2 GB | ✅ direct |
| 300-600 aa | ~3-4 GB | ✅ direct |
| >600 aa | ~5+ GB | ✅ should still fit |

ESM-IF is lightweight (~142M params). Not a VRAM concern on 12 GB.

## Phase 4: Quality Gates

### 4.1 Hard Gates (Version-Adaptive)

```python
#!/opt/conda/envs/prot/bin/python
"""Quality gates for designed sequences. Handles both NEW and LEGACY formats."""

def sequence_design_hard_gates(designs_df, wt_row, format_style):
    """Apply hard gates to LigandMPNN designs.

    Args:
        designs_df: DataFrame from mpnn_output_to_df (excluding wildtype)
        wt_row: wildtype row (Series) from the same DataFrame
        format_style: "NEW" or "LEGACY" (from mpnn_output_to_df)

    Returns:
        passed_df, failed_df, gate_report
    """
    gates = []
    failed_indices = set()

    for idx, row in designs_df.iterrows():
        reasons = []

        if format_style == "NEW":
            # NEW: overall_confidence, higher = better
            wt_conf = wt_row.get("overall_confidence", 0)
            design_conf = row.get("overall_confidence", 0)
            # Gate: design should not be drastically worse than WT
            if wt_conf > 0 and design_conf < wt_conf * 0.5:
                reasons.append(
                    f"overall_confidence {design_conf:.3f} < 0.5× WT ({wt_conf:.3f})")

        elif format_style == "LEGACY":
            # LEGACY: score, lower = better
            wt_score = wt_row.get("score", 0)
            design_score = row.get("score", 999)
            if wt_score > 0 and design_score > wt_score * 1.5:
                reasons.append(
                    f"score {design_score:.3f} > 1.5× WT ({wt_score:.3f})")

        # Recovery gate (both formats)
        rec = row.get("seq_recovery", None)
        if rec is not None:
            if rec < 0.05:
                reasons.append(f"seq_recovery {rec:.3f} < 0.05 (too divergent)")

        if reasons:
            failed_indices.add(idx)
            gates.append({"index": idx, "reasons": reasons})

    passed = designs_df[~designs_df.index.isin(failed_indices)]
    failed = designs_df[designs_df.index.isin(failed_indices)]

    report = {
        "n_input": len(designs_df),
        "n_passed": len(passed),
        "n_failed": len(failed),
        "pass_rate_pct": 100 * len(passed) / len(designs_df) if len(designs_df) > 0 else 0,
        "failures": gates,
    }

    return passed, failed, report
```

### 4.2 Soft Scores (For Ranking)

| Metric | Source | Direction | Use |
|--------|--------|-----------|-----|
| MPNN confidence/score | LigandMPNN output | NEW: higher better / LEGACY: lower better | Primary structural fit |
| ESM-IF avg_ll | ESM-IF scoring | Higher better | Cross-validation |
| Δ vs WT | design metric - WT metric | Improvement direction | Relative improvement |
| Interface recovery | Computed from sequences | Task-dependent | Diversity gauge |
| Mutation count | Sequence comparison | Lower preferred (default) | Risk minimization |

### 4.3 Cross-Validation Gate

```python
def cross_validate_mpnn_esmif(mpnn_better_than_wt, esmif_delta_ll):
    """Check if MPNN and ESM-IF agree on whether a design improves on WT.

    Args:
        mpnn_better_than_wt: bool — True if MPNN says design is better
        esmif_delta_ll: float — (design avg_ll - wt avg_ll), positive = better

    Returns:
        dict with agreement status and interpretation
    """
    esmif_says_better = esmif_delta_ll > 0.05

    if mpnn_better_than_wt and esmif_says_better:
        return {"agreement": "BOTH_BETTER", "confidence": "HIGH",
                "action": "Proceed to stability/folding evaluation"}
    elif not mpnn_better_than_wt and not esmif_says_better:
        return {"agreement": "BOTH_WORSE", "confidence": "HIGH",
                "action": "Reject — both models say worse than WT"}
    else:
        return {"agreement": "DISAGREE", "confidence": "LOW",
                "action": "Flag for manual review. Check: (1) structure quality "
                          "at mutated positions, (2) whether one model is more "
                          "reliable for this fold type, (3) proceed with caution"}
```

## Phase 5: Decision Tree — Full Workflow

```
You have a backbone structure and want to design sequences.
│
├── Step 0: Run Phase 0 environment verification
│   └── Detect CLI style (NEW/LEGACY), verify deps
│
├── Step 1: Is the structure trustworthy? (→ use prot-structure skill)
│   ├── NO → Fix structure first. Do NOT proceed.
│   └── YES ↓
│
├── Step 2: Define design scope
│   ├── Redesign everything? → no --redesigned_residues flag
│   ├── Interface only? → Extract interface residues (Phase 1.2)
│   │                     → --redesigned_residues "B17 B18 ..."
│   ├── Stability optimization? → --fixed_residues "<functional sites>"
│   └── Symmetric complex? → --tied_positions
│
├── Step 3: Choose model weights (Phase 2.4 decision tree)
│
├── Step 4: Choose temperature
│   ├── Conservative: 0.1, 50 sequences
│   ├── Balanced: 0.2, 100 sequences
│   └── Explorative: 0.3, 200 sequences
│
├── Step 5: Run LigandMPNN (Phase 2, use correct CLI style!)
│
├── Step 6: Parse output (Phase 2.6, version-adaptive parser)
│   └── CRITICAL: check score direction matches your format
│
├── Step 7: Hard gates (Phase 4.1)
│
├── Step 8: ESM-IF cross-validation on Top N (Phase 3)
│   └── Cross-validate gate (Phase 4.3)
│
├── Step 9: Rank by composite metric
│   └── Weight: 0.4 × MPNN + 0.3 × ESMIF_ll + 0.2 × mutation_count + 0.1 × recovery
│
└── Step 10: Select Top-N for downstream
    └── Folding validation: Top 20-50 / Stability: Top 10-20 / Experiment: Top 5-10
```

## Phase 6: VRAM Budget & Timing

| Task | VRAM | Wall time (<GPU_REDACTED> 12GB) | Notes |
|------|------|--------------------------|-------|
| LigandMPNN, 50 seqs, 300 aa | ~1.5 GB | ~30 sec | Very fast |
| LigandMPNN, 200 seqs, 300 aa | ~1.5 GB | ~2 min | Scales linearly |
| LigandMPNN, 50 seqs, 1000 aa | ~3 GB | ~2 min | Longer sequences use more |
| ESM-IF scoring, 1 seq, 300 aa | ~1.5 GB | ~5 sec | Per sequence |
| ESM-IF scoring, 50 seqs, 300 aa | ~1.5 GB | ~4 min | Sequential loop |

> **These are estimates.** On first real run, note actual VRAM/time and update this table.

## Failure Modes

1. **Designing on a bad backbone.** If input structure has low pLDDT loops or clashes, LigandMPNN will confidently produce sequences optimized for the wrong geometry. Always validate structure first (prot-structure skill).

2. **Forgetting to compare to wildtype.** A design score means nothing in isolation. Always parse wildtype from the first FASTA record and compute deltas.

3. **Wrong score direction.** NEW format: `overall_confidence` higher=better. LEGACY: `score` lower=better. Mixing these up silently inverts your ranking. Phase 0 detects this.

4. **Wrong chain in --redesigned_residues.** Swapping which chain you design vs fix produces plausible-looking but meaningless results. Triple-check chain IDs against the PDB.

5. **Temperature too low → near-identical sequences.** At temp=0.1 with few designable positions, expect many duplicate sequences. This is normal, not an error. If you need diversity, increase temperature.

6. **Temperature too high → sequences don't fold.** At temp≥0.4, validate with ESMFold that designs still fold to target structure.

7. **Not using SolubleMPNN when expression matters.** Standard ProteinMPNN doesn't optimize for solubility. Switch to `solublempnn_v_48_020.pt` if you've had inclusion body issues.

8. **Ignoring MPNN-ESM-IF disagreement.** Disagreement = uncertainty. Don't ignore it.

9. **Missing dependencies for ESM-IF.** ESM-IF requires PyG (torch_geometric + scatter/sparse/cluster). Phase 0 checks this. Without PyG, ESM-IF import silently fails or gives cryptic errors.

10. **biotite API mismatch.** `filter_backbone` vs `filter_peptide_backbone` — the scoring function in Phase 3 patches this automatically, but if you write your own ESM-IF code, be aware.

## One-Sentence Rule

**Design on validated structures, compare every score to wildtype, check your score direction (NEW vs LEGACY), and cross-validate LigandMPNN with ESM-IF before advancing candidates.**
