---
name: prot-structure
description: Structure prediction (ESMFold), evaluation (pLDDT/PAE/clashes/Ramachandran), and preparation for downstream protein design. Provides hard decision boundaries for when a predicted or experimental structure is trustworthy enough to use as input for sequence design, stability scoring, or interface analysis. Triggers whenever you receive a PDB, run a structure predictor, or need to decide "is this structure good enough to design on?"
homepage: https://github.com/facebookresearch/esm
metadata: { "openclaw": { "emoji": "🏗️", "requires": { "bins": ["python3"], "python": ["biopython", "numpy", "pandas", "requests", "tmtools"] } } }
---

# Protein Structure: Prediction, Evaluation & Preparation

Every downstream protein design step — sequence design, stability scoring, interface analysis — is only as good as the input structure. This skill gives you hard gates and diagnostic tools to answer: "Can I design on this structure?"

## When to Use

- You received a **PDB file** and need to assess quality before designing
- You need to **predict a structure** from sequence (ESMFold)
- You need to **validate** that a designed sequence folds to the target structure
- You need to decide whether a **predicted structure region** is trustworthy enough for design
- You need to **compare** two structures (designed vs target, ESMFold vs experimental)

## Core Philosophy

1. **Structure quality is the foundation.** A brilliant sequence design on a bad backbone is worse than no design at all — it wastes compute and gives false confidence. Gate structure quality BEFORE any downstream step.
2. **pLDDT tells you about local confidence, PAE tells you about global assembly.** They answer different questions. A protein can have high pLDDT everywhere but unreliable domain-domain orientation (high PAE). For single-domain design, pLDDT dominates. For interface/multi-domain work, PAE dominates.
3. **Always compare to a reference.** A pLDDT of 75 means nothing without context. Compare predicted structure to experimental (if available) via TM-score/RMSD. Compare designed-sequence fold-back to the target backbone.
4. **Diagnose before you design.** When structure quality is marginal, spend 5 minutes diagnosing (clashes, Ramachandran, missing residues) rather than 5 hours designing on garbage.
5. **ESMFold is fast screening, not ground truth.** Use it for rapid fold-back validation of designed sequences. For high-stakes structure assessment, cross-validate with AF2/ColabFold when possible.

## Phase 0: Environment Verification

```python
#!/opt/conda/envs/prot/bin/python
"""Phase 0: Verify structure tools are available.

ESMFold prediction strategy:
  1. DEFAULT: ESMFold API (api.esmatlas.com) — zero local deps, zero VRAM
  2. FALLBACK: Local ESMFold — requires openfold + deepspeed (hard to install)

This Phase 0 tests which prediction backend is available.
"""
import importlib, requests

print("=== prot-structure Phase 0 ===\n")

# 1. Core dependencies (evaluation tools — always needed)
deps = {
    "numpy": None, "pandas": None,
    "Bio.PDB": "biopython", "tmtools": None, "requests": None,
}
missing = []
for mod, pkg in deps.items():
    try:
        importlib.import_module(mod.split(".")[0])
    except ImportError:
        missing.append(pkg or mod)

if missing:
    print(f"❌ Missing: {missing}")
    print(f"   pip install {' '.join(missing)}")
else:
    print("✅ Core dependencies OK")

# 2. ESMFold API test (DEFAULT — preferred)
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
api_ok = False
try:
    r = requests.post(ESMFOLD_API, data="MAAAA", timeout=60)
    if r.status_code == 200 and "ATOM" in r.text:
        api_ok = True
        print("✅ ESMFold API reachable (DEFAULT prediction backend)")
    else:
        print(f"⚠️  ESMFold API returned status {r.status_code}")
except Exception as e:
    print(f"⚠️  ESMFold API unreachable: {e}")

# 3. Local ESMFold test (FALLBACK)
local_ok = False
if not api_ok:
    try:
        import esm
        model = esm.pretrained.esmfold_v1()
        local_ok = True
        print("✅ Local ESMFold available (fallback)")
        del model
    except Exception as e:
        print(f"⚠️  Local ESMFold unavailable: {e}")

if not api_ok and not local_ok:
    print("❌ NO ESMFold backend available. Cannot predict structures.")
    print("   Fix: ensure network access to api.esmatlas.com")
    print("   Or: install openfold (complex, see README)")

# 4. GPU (for evaluation tools — optional)
try:
    import torch
    if torch.cuda.is_available():
        vram = torch.cuda.get_device_properties(0).total_memory / 1024**3
        print(f"✅ GPU: {torch.cuda.get_device_name(0)} ({vram:.1f} GB)")
    else:
        print("ℹ️  No GPU (not needed for API-based prediction)")
except ImportError:
    print("ℹ️  torch not installed (not needed for API-based prediction)")

PREDICTION_BACKEND = "API" if api_ok else ("LOCAL" if local_ok else "NONE")
print(f"\n→ PREDICTION_BACKEND = {PREDICTION_BACKEND}")
print("\n=== Phase 0 complete ===")
```

## Phase 1: Structure Acquisition

### 1.1 Download Experimental Structure from PDB

```python
#!/opt/conda/envs/prot/bin/python
"""Download and inspect a PDB structure."""
import requests
from Bio.PDB import PDBParser

def download_pdb(pdb_id, out_path):
    """Download PDB from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with open(out_path, "w") as f:
        f.write(r.text)
    print(f"Downloaded {pdb_id} → {out_path}")

def inspect_pdb(pdb_path):
    """Quick inspection: chains, residues per chain, resolution if available."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    info = {"chains": {}}
    for chain in model:
        residues = [r for r in chain.get_residues() if r.id[0] == " "]
        hetatm = [r for r in chain.get_residues() if r.id[0] != " " and r.get_resname() != "HOH"]
        water = [r for r in chain.get_residues() if r.get_resname() == "HOH"]
        info["chains"][chain.id] = {
            "n_residues": len(residues),
            "n_hetatm": len(hetatm),
            "n_water": len(water),
            "first_resid": residues[0].id[1] if residues else None,
            "last_resid": residues[-1].id[1] if residues else None,
        }

    # Try to extract resolution from header
    try:
        header = structure.header
        info["resolution"] = header.get("resolution", "unknown")
        info["method"] = header.get("structure_method", "unknown")
    except:
        pass

    print(f"PDB: {pdb_path}")
    for cid, cinfo in info["chains"].items():
        print(f"  Chain {cid}: {cinfo['n_residues']} residues "
              f"({cinfo['first_resid']}-{cinfo['last_resid']}), "
              f"{cinfo['n_hetatm']} HETATM, {cinfo['n_water']} water")
    if "resolution" in info:
        print(f"  Resolution: {info['resolution']}  Method: {info['method']}")
    return info
```

### 1.2 Predict Structure with ESMFold

```python
#!/opt/conda/envs/prot/bin/python
"""Predict structure from sequence using ESMFold.

Two backends:
  1. API (default): https://api.esmatlas.com — zero VRAM, works on any machine
  2. Local (fallback): requires openfold + deepspeed (hard to install)

API limits: sequences up to ~400 aa reliably; longer may timeout.
For >400 aa, consider ColabFold or splitting into domains.
"""
import requests, os
import numpy as np
from Bio.PDB import PDBParser

ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def predict_structure_api(sequence, out_pdb, timeout=120):
    """Predict structure via ESMFold API (DEFAULT method).

    Args:
        sequence: amino acid string
        out_pdb: output PDB file path
        timeout: request timeout in seconds

    Returns:
        dict with pdb_path, plddt_mean, plddt_per_residue
    """
    print(f"Predicting via ESMFold API ({len(sequence)} aa)...")
    r = requests.post(ESMFOLD_API, data=sequence, timeout=timeout)

    if r.status_code != 200:
        raise RuntimeError(f"ESMFold API error: HTTP {r.status_code}\n{r.text[:200]}")

    pdb_text = r.text
    if "ATOM" not in pdb_text:
        raise RuntimeError(f"ESMFold API returned invalid PDB:\n{pdb_text[:200]}")

    with open(out_pdb, "w") as f:
        f.write(pdb_text)

    plddt = extract_plddt_from_pdb(out_pdb)

    result = {
        "pdb_path": out_pdb,
        "sequence_length": len(sequence),
        "plddt_mean": float(plddt.mean()) if len(plddt) > 0 else 0,
        "plddt_min": float(plddt.min()) if len(plddt) > 0 else 0,
        "plddt_per_residue": plddt,
        "backend": "API",
    }

    print(f"  → {out_pdb} | Mean pLDDT: {result['plddt_mean']:.1f} | "
          f"Min pLDDT: {result['plddt_min']:.1f}")
    return result


def predict_structure_local(sequence, out_pdb, device="cuda", max_len_fp32=300):
    """Predict structure via local ESMFold (FALLBACK — requires openfold).

    Only use if API is unavailable. Requires:
    - pip install fair-esm[esmfold] openfold deepspeed
    - Matching versions (torch/openfold/deepspeed compatibility is fragile)

    VRAM budget (<GPU_REDACTED> 12GB):
      ≤300 aa:  ~6-8 GB  → fp32
      300-400:  ~8-10 GB → fp16
      >400:     >10 GB   → CPU or API
    """
    import torch, esm

    model = esm.pretrained.esmfold_v1()
    model = model.eval().to(device)

    use_fp16 = len(sequence) > max_len_fp32 and device == "cuda"

    with torch.no_grad():
        if use_fp16:
            with torch.cuda.amp.autocast():
                output = model.infer_pdb(sequence)
        else:
            output = model.infer_pdb(sequence)

    with open(out_pdb, "w") as f:
        f.write(output)

    plddt = extract_plddt_from_pdb(out_pdb)
    del model
    torch.cuda.empty_cache()

    return {
        "pdb_path": out_pdb,
        "sequence_length": len(sequence),
        "plddt_mean": float(plddt.mean()) if len(plddt) > 0 else 0,
        "plddt_min": float(plddt.min()) if len(plddt) > 0 else 0,
        "plddt_per_residue": plddt,
        "backend": "LOCAL",
    }


def predict_structure(sequence, out_pdb, backend="API", **kwargs):
    """Predict structure — auto-selects backend.

    Args:
        sequence: amino acid string
        out_pdb: output PDB path
        backend: "API" (default), "LOCAL", or "AUTO"

    Returns:
        dict with prediction results
    """
    if backend == "AUTO":
        try:
            return predict_structure_api(sequence, out_pdb)
        except Exception as e:
            print(f"⚠️  API failed ({e}), trying local...")
            return predict_structure_local(sequence, out_pdb, **kwargs)
    elif backend == "API":
        return predict_structure_api(sequence, out_pdb)
    elif backend == "LOCAL":
        return predict_structure_local(sequence, out_pdb, **kwargs)
    else:
        raise ValueError(f"Unknown backend: {backend}")


def extract_plddt_from_pdb(pdb_path):
    """Extract per-residue pLDDT from B-factor column of ESMFold output.

    Handles unit mismatch: ESMFold API returns B-factors in 0-1 range,
    while local ESMFold uses 0-100. This function auto-detects and
    normalizes to 0-100 scale.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    plddt = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            ca = residue["CA"] if "CA" in residue else None
            if ca is not None:
                plddt.append(ca.get_bfactor())

    plddt = np.array(plddt)

    # Auto-detect unit: if max ≤ 1.0, values are 0-1 → scale to 0-100
    if len(plddt) > 0 and plddt.max() <= 1.0:
        plddt = plddt * 100.0

    return plddt
```

### 1.3 Prediction Backend Strategy

| Backend | Pros | Cons | When to use |
|---------|------|------|-------------|
| **API** (default) | Zero VRAM, zero deps, fast | Needs network; ~400 aa limit; rate limits | Always try first |
| **Local** (fallback) | No network needed, no length limit | openfold/deepspeed install is fragile; 12GB VRAM limit | Only if API down |

**API rate limits**: If you get HTTP 429 or timeouts, add `time.sleep(2)` between requests.

**Sequence length**: API handles up to ~400 aa reliably. For >400 aa, split into domains or use ColabFold.

```python
def predict_structure_safe(sequence, out_pdb, backend="API"):
    """VRAM-safe + network-safe wrapper. Tries API first, falls back to local."""
    if len(sequence) > 400 and backend == "API":
        print(f"⚠️  Sequence {len(sequence)} aa — API may timeout. Trying anyway...")
    return predict_structure(sequence, out_pdb, backend="AUTO")
```

## Phase 2: Structure Evaluation — Decision Boundaries

### 2.1 The Decision Tree (Hard Gates)

```
You have a structure (experimental or predicted). Can you design on it?
│
├── Step 1: Global pLDDT check (predicted structures only)
│   ├── Mean pLDDT ≥ 85 → HIGH confidence. Proceed.
│   ├── Mean pLDDT 75-85 → MODERATE. Proceed with local checks (Step 2).
│   ├── Mean pLDDT 60-75 → LOW. Only proceed if:
│   │   - Design region itself has pLDDT ≥ 75, AND
│   │   - You freeze all low-confidence regions
│   └── Mean pLDDT < 60 → REJECT. Do not design. Structure is unreliable.
│
├── Step 2: Local checks on design region
│   ├── All residues in design region pLDDT ≥ 70? → OK
│   ├── Contiguous stretch ≥10 residues with pLDDT < 70?
│   │   ├── Region is non-functional surface loop → FREEZE it, design around it
│   │   ├── Region is at interface/functional site → DO NOT design.
│   │   │   Try: (a) AF2/ColabFold rebuild, (b) template-based modeling
│   │   └── Region is in core → REJECT structure for design
│   └── Scattered 1-3 residues pLDDT < 70 → OK, but freeze those positions
│
├── Step 3: Clash and geometry check (all structures)
│   ├── Severe clashes (>10 atom pairs < 1.5 Å) → REJECT or re-minimize
│   ├── Ramachandran outliers > 5% → WARNING (check if in design region)
│   └── Missing heavy atoms in design region → REJECT (incomplete structure)
│
├── Step 4: Multi-chain / interface check (complexes only)
│   ├── Inter-chain PAE < 10 Å (or ipTM > 0.7 for AF2-multimer) → Interface OK
│   ├── Inter-chain PAE 10-20 Å → Interface UNCERTAIN. Design with caution.
│   └── Inter-chain PAE > 20 Å (or ipTM < 0.5) → Interface UNRELIABLE.
│       Do NOT do interface design. Seek experimental complex structure.
│
└── Step 5: Sanity checks
    ├── Repeat prediction 2-3x: backbone RMSD < 1 Å → Stable prediction
    ├── Secondary structure makes biological sense (not all-coil for a known β-barrel)
    └── Hydrophobic core exists (not inside-out)
```

### 2.2 Evaluation Code

```python
#!/opt/conda/envs/prot/bin/python
"""Structure evaluation: pLDDT analysis, clash detection, Ramachandran."""
import numpy as np
from Bio.PDB import PDBParser

def evaluate_plddt(pdb_path, design_residues=None):
    """Evaluate pLDDT distribution and apply decision gates.

    Args:
        pdb_path: PDB with pLDDT in B-factor (ESMFold output)
        design_residues: list of (chain, resseq) tuples for design region
                         If None, evaluates globally.

    Returns:
        dict with assessment and per-region stats
    """
    plddt = extract_plddt_from_pdb(pdb_path)

    # Global stats
    result = {
        "global_mean": float(plddt.mean()),
        "global_min": float(plddt.min()),
        "global_median": float(np.median(plddt)),
        "pct_above_90": float((plddt >= 90).mean() * 100),
        "pct_above_70": float((plddt >= 70).mean() * 100),
        "pct_below_60": float((plddt < 60).mean() * 100),
    }

    # Global gate
    mean = result["global_mean"]
    if mean >= 85:
        result["global_gate"] = "HIGH"
        result["global_action"] = "Proceed to design"
    elif mean >= 75:
        result["global_gate"] = "MODERATE"
        result["global_action"] = "Proceed with local checks on design region"
    elif mean >= 60:
        result["global_gate"] = "LOW"
        result["global_action"] = "Only proceed if design region pLDDT ≥ 75; freeze low regions"
    else:
        result["global_gate"] = "REJECT"
        result["global_action"] = "Do NOT design. Structure unreliable."

    # Design region analysis (if specified)
    if design_residues is not None:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("s", pdb_path)
        model = list(structure.get_models())[0]

        design_plddt = []
        design_set = set(design_residues)
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                if (chain.id, res.id[1]) in design_set:
                    ca = res["CA"] if "CA" in res else None
                    if ca:
                        design_plddt.append(ca.get_bfactor())

        if design_plddt:
            design_plddt = np.array(design_plddt)
            result["design_region_mean"] = float(design_plddt.mean())
            result["design_region_min"] = float(design_plddt.min())
            result["design_region_pct_below_70"] = float((design_plddt < 70).mean() * 100)

            # Find contiguous low-confidence stretches
            low_mask = design_plddt < 70
            max_contiguous_low = _max_contiguous_true(low_mask)
            result["design_region_max_contiguous_low"] = max_contiguous_low

            # Design region gate
            if design_plddt.mean() >= 75 and max_contiguous_low < 10:
                result["design_gate"] = "OK"
                result["design_action"] = "Design region quality sufficient"
            elif max_contiguous_low >= 10:
                result["design_gate"] = "REJECT_REGION"
                result["design_action"] = (
                    f"Contiguous low-confidence stretch ({max_contiguous_low} residues) "
                    "in design region. Freeze region or rebuild with AF2."
                )
            else:
                result["design_gate"] = "CAUTION"
                result["design_action"] = "Freeze low-pLDDT positions; design others"

    return result

def _max_contiguous_true(mask):
    """Find longest contiguous True stretch in boolean array."""
    max_len = 0
    current = 0
    for v in mask:
        if v:
            current += 1
            max_len = max(max_len, current)
        else:
            current = 0
    return max_len


def detect_clashes(pdb_path, clash_cutoff=1.5, min_seq_sep=3):
    """Detect steric clashes (heavy atom pairs closer than cutoff).

    Properly excludes bonded neighbors: atoms in residues within
    min_seq_sep sequence positions are skipped (they have legitimate
    close contacts from covalent bonds and angles).

    Args:
        pdb_path: PDB file
        clash_cutoff: distance in Å below which atoms are clashing (default 1.5)
        min_seq_sep: minimum residue sequence separation to consider (default 3,
                     skips i/i+1/i+2 which are bonded or angle-connected)

    Returns:
        dict with clash count, assessment, worst pairs, and peptide flag
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    atoms = []
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            resseq = res.id[1]
            for atom in res:
                if atom.element == "H":
                    continue
                atoms.append({
                    "coord": atom.coord,
                    "name": atom.get_name(),
                    "chain": chain.id,
                    "resseq": resseq,
                    "resname": res.get_resname(),
                    "label": f"{chain.id}{resseq}{res.get_resname()}",
                })

    n_residues = len(set((a["chain"], a["resseq"]) for a in atoms))
    is_short_peptide = n_residues < 30

    # Check all pairs with sufficient sequence separation
    clashes = []
    coords = np.array([a["coord"] for a in atoms])
    n = len(coords)

    for i in range(n):
        for j in range(i + 1, n):
            # Skip same residue
            if atoms[i]["label"] == atoms[j]["label"]:
                continue
            # Skip bonded neighbors (same chain, close in sequence)
            if atoms[i]["chain"] == atoms[j]["chain"]:
                seq_sep = abs(atoms[i]["resseq"] - atoms[j]["resseq"])
                if seq_sep < min_seq_sep:
                    continue
            d = np.linalg.norm(coords[i] - coords[j])
            if d < clash_cutoff:
                clashes.append({
                    "atom1": f"{atoms[i]['label']}.{atoms[i]['name']}",
                    "atom2": f"{atoms[j]['label']}.{atoms[j]['name']}",
                    "distance": float(d),
                })

    result = {
        "n_clashes": len(clashes),
        "clash_cutoff": clash_cutoff,
        "n_residues": n_residues,
        "is_short_peptide": is_short_peptide,
    }

    # Scale gate thresholds for short peptides
    severe_threshold = 3 if is_short_peptide else 10

    if len(clashes) > severe_threshold:
        result["gate"] = "REJECT"
        result["action"] = (f"{len(clashes)} clashes (>{severe_threshold} for "
                           f"{'peptide' if is_short_peptide else 'protein'}). "
                           f"Re-minimize or reject.")
    elif len(clashes) > 0:
        result["gate"] = "WARNING"
        result["action"] = f"{len(clashes)} minor clashes. Check if in design region."
    else:
        result["gate"] = "OK"
        result["action"] = "No clashes detected."

    result["worst_clashes"] = sorted(clashes, key=lambda x: x["distance"])[:5]

    return result
```

## Phase 3: Structure Comparison (Design Validation)

### 3.1 Compare Two Structures via TM-score

```python
#!/opt/conda/envs/prot/bin/python
"""Compare two structures using TM-score and RMSD."""
import tmtools
import numpy as np
from Bio.PDB import PDBParser

def extract_ca_coords(pdb_path, chain_id=None):
    """Extract CA coordinates and sequence from a PDB.

    Args:
        pdb_path: PDB file path
        chain_id: specific chain (None = first chain)

    Returns:
        coords (N×3 numpy array), sequence (str)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    if chain_id:
        chains = [model[chain_id]]
    else:
        chains = [list(model.get_chains())[0]]

    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }

    coords = []
    seq = []
    for chain in chains:
        for res in chain:
            if res.id[0] != " ":
                continue
            if "CA" not in res:
                continue
            coords.append(res["CA"].coord)
            seq.append(three_to_one.get(res.get_resname(), "X"))

    return np.array(coords), "".join(seq)


def compare_structures(pdb_a, pdb_b, chain_a=None, chain_b=None):
    """Compare two structures using TM-score and CA RMSD.

    Use cases:
    - Design validation: compare ESMFold(designed_seq) vs target backbone
    - Prediction quality: compare ESMFold vs experimental structure

    Returns:
        dict with tm_score, rmsd, aligned_length, and assessment
    """
    coords_a, seq_a = extract_ca_coords(pdb_a, chain_a)
    coords_b, seq_b = extract_ca_coords(pdb_b, chain_b)

    # TM-align
    result = tmtools.tm_align(coords_a, coords_b, seq_a, seq_b)

    tm_score = result.tm_norm_chain1  # normalized by length of chain 1
    rmsd = result.rmsd

    # Assessment
    if tm_score >= 0.9:
        assessment = "EXCELLENT"
        action = "Structures are highly similar. Design likely folds correctly."
    elif tm_score >= 0.7:
        assessment = "GOOD"
        action = "Structures are similar. Minor deviations acceptable."
    elif tm_score >= 0.5:
        assessment = "MARGINAL"
        action = "Significant deviations. Check which regions differ."
    else:
        assessment = "POOR"
        action = "Structures are dissimilar. Design likely does NOT fold to target."

    comparison = {
        "tm_score": float(tm_score),
        "rmsd": float(rmsd),
        "aligned_length": int(result.aligned_length) if hasattr(result, 'aligned_length') else min(len(seq_a), len(seq_b)),
        "len_a": len(seq_a),
        "len_b": len(seq_b),
        "assessment": assessment,
        "action": action,
    }

    print(f"TM-score: {tm_score:.3f} | RMSD: {rmsd:.2f} Å | Assessment: {assessment}")

    return comparison
```

### 3.2 Wildtype Sanity Check (Run BEFORE Design Fold-Back)

```python
#!/opt/conda/envs/prot/bin/python
"""Check if ESMFold can reliably predict the wildtype structure.

If ESMFold cannot fold the WILDTYPE well, fold-back results for
designs are meaningless. This check must run BEFORE batch fold-back.

Discovered in IL-2 (1Z92): ESMFold pLDDT≈45, TM≈0.49 for WT.
Proteins with disulfide bonds, multi-domain architectures, or
unusual folds may score poorly. In those cases, fold-back should
be downgraded from hard gate to diagnostic-only.
"""
import os

def wildtype_sanity_check(wt_sequence, target_pdb, chain_id=None,
                           work_dir="foldback_wt_sanity", backend="API"):
    """Predict wildtype structure and compare to experimental.

    Args:
        wt_sequence: wildtype amino acid string
        target_pdb: experimental structure PDB
        chain_id: chain to compare
        work_dir: temp directory
        backend: "API" (default) or "LOCAL"

    Returns:
        dict with WT fold-back results + recommended gate_mode
    """
    os.makedirs(work_dir, exist_ok=True)

    # Predict WT
    wt_pdb = os.path.join(work_dir, "wt_esmfold.pdb")
    wt_pred = predict_structure(wt_sequence, wt_pdb, backend=backend)

    # Compare to experimental
    comparison = compare_structures(target_pdb, wt_pdb, chain_a=chain_id)

    wt_tm = comparison["tm_score"]
    wt_plddt = wt_pred["plddt_mean"]

    # Determine gate mode
    if wt_plddt >= 60 and wt_tm >= 0.7:
        gate_mode = "HARD_GATE"
        recommendation = ("ESMFold folds WT reliably. "
                         "Fold-back can be used as a hard gate.")
    elif wt_plddt >= 45 and wt_tm >= 0.4:
        gate_mode = "DIAGNOSTIC"
        recommendation = (f"ESMFold partially folds WT (pLDDT={wt_plddt:.1f}, "
                         f"TM={wt_tm:.2f}). Fold-back downgraded to diagnostic: "
                         f"report but do not hard-reject. Compare designs RELATIVE "
                         f"to WT fold-back, not against absolute thresholds.")
    else:
        gate_mode = "SKIP"
        recommendation = (f"ESMFold cannot fold WT (pLDDT={wt_plddt:.1f}, "
                         f"TM={wt_tm:.2f}). Fold-back is unreliable for this "
                         f"target. Skip fold-back gate entirely.")

    result = {
        "wt_plddt": wt_plddt,
        "wt_tm": wt_tm,
        "wt_rmsd": comparison["rmsd"],
        "gate_mode": gate_mode,
        "recommendation": recommendation,
        "wt_pred_pdb": wt_pdb,
    }

    print(f"\n=== WT Sanity Check ===")
    print(f"  WT pLDDT: {wt_plddt:.1f}")
    print(f"  WT TM-score: {wt_tm:.3f}")
    print(f"  Gate mode: {gate_mode}")
    print(f"  {recommendation}")

    return result
```

### 3.2.1 Gate Mode Thresholds

| WT pLDDT | WT TM-score | Gate Mode | Fold-back behavior |
|----------|-------------|-----------|-------------------|
| ≥ 60 | ≥ 0.7 | **HARD_GATE** | Normal: PASS/MARGINAL/FAIL as usual |
| 45-60 | 0.4-0.7 | **DIAGNOSTIC** | Report only; compare design vs WT (relative), not absolute |
| < 45 | < 0.4 | **SKIP** | Do not use fold-back for this target |

> **Why?** ESMFold has known weaknesses: disulfide-rich proteins (IL-2), multi-domain proteins, membrane proteins, IDPs. If it can't fold the wildtype, it can't validate designs. The WT sanity check prevents false rejections like the IL-2 1Z92 case where ALL designs (including correct ones) would be rejected because ESMFold simply can't predict this fold well.

### 3.3 Fold-Back Validation (Design → ESMFold → Compare)

```python
#!/opt/conda/envs/prot/bin/python
"""Validate that a designed sequence folds to the target structure."""
import os

def foldback_validation(designed_sequence, target_pdb, chain_id=None,
                         work_dir="foldback_tmp", backend="API",
                         gate_mode="HARD_GATE", wt_sanity=None):
    """Predict structure of designed sequence and compare to target.

    This is THE critical quality check for sequence design:
    if ESMFold(designed_seq) ≠ target backbone, the design is likely wrong.

    Handles three modes:
    - HARD_GATE: standard pass/fail (default, when WT sanity passes)
    - DIAGNOSTIC: report metrics but don't reject (when WT sanity is marginal)
    - SKIP: skip entirely (when ESMFold can't fold the WT)

    Handles short peptides (<30 aa) differently:
    - TM-score is statistically unreliable for short sequences
    - Uses CA RMSD as primary metric instead
    - pLDDT gate still applies

    Args:
        designed_sequence: amino acid string
        target_pdb: path to target backbone PDB
        chain_id: chain to compare in target
        work_dir: temp directory for predicted PDB
        backend: "API" (default) or "LOCAL"
        gate_mode: "HARD_GATE", "DIAGNOSTIC", or "SKIP" (from wt_sanity_check)
        wt_sanity: result from wildtype_sanity_check() (for relative comparison)

    Returns:
        dict with prediction result + comparison + pass/fail gate
    """
    if gate_mode == "SKIP":
        return {
            "prediction": None,
            "comparison": None,
            "gate": "SKIPPED",
            "action": "Fold-back skipped: ESMFold cannot reliably fold this target.",
            "gate_mode": gate_mode,
        }

    os.makedirs(work_dir, exist_ok=True)

    # Predict
    pred_pdb = os.path.join(work_dir, "designed_fold.pdb")
    pred_result = predict_structure(designed_sequence, pred_pdb, backend=backend)

    # Compare
    comparison = compare_structures(target_pdb, pred_pdb, chain_a=chain_id)

    # Gate — different thresholds for short peptides vs proteins
    tm = comparison["tm_score"]
    rmsd = comparison["rmsd"]
    plddt = pred_result["plddt_mean"]
    is_short = len(designed_sequence) < 30

    if gate_mode == "DIAGNOSTIC":
        # Diagnostic mode: compare relative to WT, not absolute thresholds
        gate = "DIAGNOSTIC"
        if wt_sanity:
            delta_tm = tm - wt_sanity["wt_tm"]
            delta_plddt = plddt - wt_sanity["wt_plddt"]
            action = (f"DIAGNOSTIC (WT unreliable): design TM={tm:.2f} "
                     f"(ΔWT={delta_tm:+.2f}), pLDDT={plddt:.1f} "
                     f"(ΔWT={delta_plddt:+.1f}). "
                     f"{'Design ≥ WT' if delta_tm >= -0.05 and delta_plddt >= -5 else 'Design < WT — flag for review'}.")
        else:
            action = (f"DIAGNOSTIC: TM={tm:.2f}, pLDDT={plddt:.1f}. "
                     f"Not gated (WT sanity check did not pass).")
    elif is_short:
        # Short peptides: TM-score unreliable, use RMSD + pLDDT
        if rmsd <= 2.0 and plddt >= 70:
            gate = "PASS"
            action = (f"Short peptide fold-back OK (RMSD={rmsd:.2f}Å, "
                      f"pLDDT={plddt:.1f}). Note: TM-score={tm:.2f} not "
                      f"reliable for {len(designed_sequence)} aa.")
        elif rmsd <= 3.5 and plddt >= 60:
            gate = "MARGINAL"
            action = (f"Short peptide partially matches (RMSD={rmsd:.2f}Å). "
                      f"Check if key secondary structure is preserved.")
        else:
            gate = "FAIL"
            action = (f"Short peptide does NOT match target "
                      f"(RMSD={rmsd:.2f}Å, pLDDT={plddt:.1f}).")
    else:
        # Standard proteins: TM-score + pLDDT
        if tm >= 0.85 and plddt >= 75:
            gate = "PASS"
            action = "Design folds to target structure with high confidence."
        elif tm >= 0.70 and plddt >= 70:
            gate = "MARGINAL"
            action = "Design partially matches. Check deviating regions."
        else:
            gate = "FAIL"
            action = (f"Design does NOT fold to target "
                      f"(TM={tm:.2f}, pLDDT={plddt:.1f}). Reject or redesign.")

    return {
        "prediction": pred_result,
        "comparison": comparison,
        "gate": gate,
        "action": action,
        "is_short_peptide": is_short,
        "primary_metric": "RMSD" if is_short else "TM-score",
        "gate_mode": gate_mode,
    }

# Usage:
# # Step 1: WT sanity check
# wt_sanity = wildtype_sanity_check("APTSSSTKKT...", "1Z92_clean.pdb", chain_id="A")
# # Step 2: Design fold-back with appropriate gate mode
# result = foldback_validation("APTSSSTKKY...", "1Z92_clean.pdb", chain_id="A",
#                               gate_mode=wt_sanity["gate_mode"], wt_sanity=wt_sanity)
```

### 3.4 Fold-Back Gate Thresholds

**Standard proteins (≥30 aa):**

| Metric | PASS | MARGINAL | FAIL |
|--------|------|----------|------|
| TM-score (designed vs target) | ≥ 0.85 | 0.70-0.85 | < 0.70 |
| Mean pLDDT (designed fold) | ≥ 75 | 70-75 | < 70 |
| Both required | TM≥0.85 AND pLDDT≥75 | TM≥0.70 AND pLDDT≥70 | either fails |

**Short peptides (<30 aa):**

| Metric | PASS | MARGINAL | FAIL |
|--------|------|----------|------|
| CA RMSD (designed vs target) | ≤ 2.0 Å | 2.0-3.5 Å | > 3.5 Å |
| Mean pLDDT (designed fold) | ≥ 70 | 60-70 | < 60 |
| Primary metric | RMSD (TM-score reported but not gated) | RMSD | RMSD |

> **Why different for short peptides?** TM-score normalizes by sequence length. For sequences <30 aa, the normalization makes TM-score statistically unstable — a 13-aa peptide can have TM=0.4 even when the RMSD is <1 Å. Use RMSD as the primary metric for short peptides.

## Phase 4: Batch Evaluation Pipeline

### 4.1 Evaluate Multiple Designs

```python
#!/opt/conda/envs/prot/bin/python
"""Batch fold-back validation for multiple designed sequences."""
import pandas as pd
import os, gc, torch

def batch_foldback(sequences, names, target_pdb, chain_id=None,
                   work_dir="foldback_batch", device="cuda"):
    """Fold-back validate a batch of designed sequences.

    Handles VRAM by processing one at a time and clearing cache.

    Args:
        sequences: list of amino acid strings
        names: list of names/labels for each sequence
        target_pdb: target backbone PDB
        chain_id: chain to compare in target
        work_dir: working directory

    Returns:
        DataFrame with results + summary
    """
    os.makedirs(work_dir, exist_ok=True)
    results = []

    for i, (seq, name) in enumerate(zip(sequences, names)):
        print(f"\n[{i+1}/{len(sequences)}] Fold-back: {name} ({len(seq)} aa)")
        try:
            pred_pdb = os.path.join(work_dir, f"{name}.pdb")
            pred = predict_structure_safe(seq, pred_pdb, device=device)
            comp = compare_structures(target_pdb, pred_pdb, chain_a=chain_id)

            tm = comp["tm_score"]
            plddt = pred["plddt_mean"]
            gate = "PASS" if (tm >= 0.85 and plddt >= 75) else \
                   "MARGINAL" if (tm >= 0.70 and plddt >= 70) else "FAIL"

            results.append({
                "name": name,
                "sequence": seq,
                "plddt_mean": plddt,
                "plddt_min": pred["plddt_min"],
                "tm_score": tm,
                "rmsd": comp["rmsd"],
                "gate": gate,
            })
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({
                "name": name, "sequence": seq,
                "plddt_mean": 0, "plddt_min": 0,
                "tm_score": 0, "rmsd": 999, "gate": "ERROR",
            })

        # Clear VRAM between predictions
        if device == "cuda":
            torch.cuda.empty_cache()
            gc.collect()

    df = pd.DataFrame(results)

    # Summary
    n_pass = (df["gate"] == "PASS").sum()
    n_marginal = (df["gate"] == "MARGINAL").sum()
    n_fail = (df["gate"] == "FAIL").sum()
    n_error = (df["gate"] == "ERROR").sum()

    print(f"\n=== Batch fold-back summary ===")
    print(f"PASS: {n_pass}  MARGINAL: {n_marginal}  FAIL: {n_fail}  ERROR: {n_error}")
    print(f"Mean TM-score: {df[df['gate']!='ERROR']['tm_score'].mean():.3f}")

    return df
```

## Phase 5: Decision Tree — Full Workflow

```
You need to evaluate or prepare a structure for protein design.
│
├── Step 1: Where does the structure come from?
│   ├── Experimental PDB → Download (Phase 1.1), inspect, go to Step 2
│   ├── ESMFold prediction → Predict (Phase 1.2), extract pLDDT, go to Step 2
│   └── Other predictor (AF2/ColabFold) → Load, extract confidence, go to Step 2
│
├── Step 2: Evaluate structure quality (Phase 2)
│   ├── Run evaluate_plddt() with design_residues
│   ├── Run detect_clashes()
│   ├── Check: do gates pass? (See decision tree in Phase 2.1)
│   ├── REJECT → Stop. Do not proceed to design.
│   ├── CAUTION → Freeze problematic regions, proceed with reduced scope
│   └── OK → Proceed to design
│
├── Step 3: (If complex) Evaluate interface reliability
│   └── Check inter-chain PAE or ipTM (Phase 2.1 Step 4)
│
├── Step 4: Proceed to downstream skill
│   ├── Sequence design → prot-seqdesign
│   ├── Stability scoring → prot-stability
│   └── Interface analysis → prot-interface
│
└── Step 5: After design — fold-back validation (Phase 3.2)
    ├── For Top-N designed sequences
    ├── Run foldback_validation() or batch_foldback()
    └── Only advance designs that PASS fold-back gate
```

## Phase 6: VRAM Budget & Timing

| Task | VRAM | Wall time (<GPU_REDACTED>) | Notes |
|------|------|---------------------|-------|
| ESMFold API, ≤300 aa | 0 GB | ~5-15 sec | Network dependent |
| ESMFold API, 300-400 aa | 0 GB | ~10-30 sec | May timeout >400 aa |
| ESMFold API, 20 designs batch | 0 GB | ~3-8 min | Add sleep(2) between requests |
| ESMFold local, ≤300 aa | ~6-8 GB | ~10-20 sec | Requires openfold |
| pLDDT evaluation | 0 GB | ~1 sec | CPU only, from PDB B-factors |
| TM-score comparison | 0 GB | ~1 sec | CPU only |
| Clash detection | 0 GB | ~2-5 sec | CPU only |

> **First real run: note actual VRAM/time and update this table.**

## Failure Modes

1. **Designing on a low-pLDDT backbone.** The #1 failure mode. LigandMPNN will happily design sequences for a backbone that doesn't exist. The structure quality gate (Phase 2) exists to prevent this.

2. **Ignoring PAE for multi-chain work.** A complex can have perfect per-residue pLDDT but unreliable domain/chain arrangement. For interface design, PAE/ipTM matters more than pLDDT.

3. **ESMFold hallucination — confident but wrong.** ESMFold can produce high-pLDDT structures for sequences that don't actually fold that way. Cross-validate by repeating prediction (different seeds if possible) and checking biological plausibility (secondary structure, hydrophobic core).

4. **Skipping fold-back validation.** After sequence design, always predict the structure of the designed sequence and compare to the target. This catches designs that score well on paper but don't actually fold to the intended structure.

5. **Using TM-score gates for short peptides (<30 aa).** TM-score is statistically unreliable for short sequences due to length normalization. A 13-aa peptide can score TM=0.4 even with RMSD <1 Å. Use CA RMSD as primary metric for peptides <30 aa (see Phase 3.3).

6. **pLDDT unit mismatch between API and local ESMFold.** The ESMFold API returns B-factors in 0-1 range; local ESMFold uses 0-100. The `extract_plddt_from_pdb()` function auto-detects and normalizes, but if you write custom code, check `max(bfactors)` to determine the scale.

7. **VRAM OOM on long sequences (local only).** ESMFold uses quadratic memory in sequence length. The API backend avoids this entirely.

8. **Comparing structures of different lengths.** TM-score handles this gracefully, but RMSD doesn't. Always use TM-score as the primary comparison metric for proteins ≥30 aa, not raw RMSD.

9. **Trusting ESMFold as ground truth.** ESMFold is fast and good but not AF2-level accurate. For critical decisions (advancing to experiments), cross-validate with AF2/ColabFold if possible.

10. **Skipping the WT sanity check.** The #1 cause of false rejections. If ESMFold can't fold the wildtype (pLDDT<60 or TM<0.7 vs experimental), ALL designs will also fail — not because they're bad, but because ESMFold can't predict this fold. Always run `wildtype_sanity_check()` before batch fold-back. Discovered in IL-2 (1Z92): ESMFold pLDDT≈45 for WT, causing every design to be falsely rejected. After adding WT sanity check, fold-back correctly downgraded to diagnostic mode.

11. **Clash detection false positives on short peptides.** The clash detector uses min_seq_sep=3 to exclude bonded neighbors, but very short peptides have fewer long-range contacts, making the absolute clash count misleading. Gate thresholds are automatically adjusted for peptides <30 aa.

## One-Sentence Rule

**Gate structure quality with hard pLDDT/PAE thresholds before any design step, always run a WT sanity check before fold-back validation, and if ESMFold can't fold the wildtype — downgrade fold-back to diagnostic rather than falsely rejecting every design.**
