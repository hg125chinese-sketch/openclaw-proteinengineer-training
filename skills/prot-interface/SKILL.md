---
name: prot-interface
description: Protein-protein interface analysis, hotspot identification, contact mapping, and binding energy estimation. Covers interface residue extraction, contact classification, per-residue energy decomposition via ESM-2, alanine scanning proxy, and interface-aware design constraints. Triggers whenever you need to analyze a protein-protein or protein-peptide binding interface, identify hotspots, or prepare constraints for interface redesign.
homepage: https://github.com/facebookresearch/esm
metadata: { "openclaw": { "emoji": "🤝", "requires": { "bins": ["python3"], "python": ["biopython", "numpy", "pandas", "requests", "torch", "fair-esm"] } } }
---

# Protein Interface: Analysis, Hotspots & Design Constraints

Interface analysis bridges structure and function. Before designing at an interface, you need to know: which residues make contact, which contacts are energetically important (hotspots), and which positions are safe to mutate. This skill gives you the tools to answer those questions and feed the results into prot-seqdesign.

## When to Use

- Analyzing a **protein-protein or protein-peptide complex** before interface design
- Identifying **hotspot residues** that drive binding
- Preparing **design constraints** for prot-seqdesign (which positions to fix vs redesign)
- Evaluating whether a **designed mutation** might disrupt binding
- Comparing **interface properties** between wildtype and designed complexes
- Deciding whether to do **interface redesign vs scaffold redesign**

## Core Philosophy

1. **Map the interface before you design it.** Blind redesign of an interface wastes compute and risks disrupting critical contacts. Spend 10 minutes mapping contacts and hotspots to save hours of failed designs.
2. **Hotspots are rare and precious.** Typically 5-10% of interface residues contribute >70% of binding energy (Bogan & Thorn, 1998). Identify them and protect them — or deliberately target them if you want to disrupt binding.
3. **Contact type matters.** A hydrogen bond across the interface is worth more than a van der Waals contact. Classify contacts by type (H-bond, salt bridge, hydrophobic, aromatic) to understand what drives binding.
4. **ESM-2 alanine scanning is a fast proxy, not ΔΔG.** It tells you how "unusual" a position is in evolutionary context — positions where alanine is very unlikely are probably functionally important. But it cannot distinguish binding-critical from fold-critical importance.
5. **Interface design is constrained optimization.** You want to improve affinity (or change specificity) while maintaining fold stability. Always check both prot-stability AND prot-interface metrics for designed mutations.

## Phase 0: Environment Verification

```python
#!/opt/conda/envs/prot/bin/python
"""Phase 0: Verify interface analysis tools."""
import importlib

print("=== prot-interface Phase 0 ===\n")

deps = {
    "Bio.PDB": "biopython",
    "numpy": None, "pandas": None,
    "torch": None, "esm": "fair-esm",
    "requests": None,
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
    print("✅ All dependencies OK")

# GPU (for ESM-2 hotspot analysis)
import torch
if torch.cuda.is_available():
    vram = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f"✅ GPU: {torch.cuda.get_device_name(0)} ({vram:.1f} GB)")
else:
    print("⚠️  No GPU — ESM-2 alanine scanning will be slow")

print("\n=== Phase 0 complete ===")
```

## Phase 1: Interface Extraction & Contact Mapping

### 1.1 Extract Interface Residues

```python
#!/opt/conda/envs/prot/bin/python
"""Extract interface residues between two chains."""
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

def extract_interface(pdb_path, chain_a, chain_b, cutoff=5.0):
    """Find residues at the interface between two chains.

    A residue is "at the interface" if any of its heavy atoms is
    within cutoff Å of any heavy atom in the other chain.

    Args:
        pdb_path: PDB file path
        chain_a: chain ID of receptor/partner A
        chain_b: chain ID of ligand/partner B
        cutoff: distance cutoff in Å (default 5.0)

    Returns:
        dict with interface residues for each chain, counts, and details
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    # Get atoms for each chain
    atoms_a = [atom for atom in model[chain_a].get_atoms() if atom.element != "H"]
    atoms_b = [atom for atom in model[chain_b].get_atoms() if atom.element != "H"]

    # Build neighbor search on chain B
    ns_b = NeighborSearch(atoms_b)

    # Find chain A residues with contacts to chain B
    interface_a = set()
    contacts = []
    for atom_a in atoms_a:
        neighbors = ns_b.search(atom_a.coord, cutoff, "A")  # atom level
        for atom_b in neighbors:
            res_a = atom_a.get_parent()
            res_b = atom_b.get_parent()
            if res_a.id[0] != " " or res_b.id[0] != " ":
                continue  # skip HETATM
            interface_a.add((chain_a, res_a.id[1], res_a.get_resname()))
            d = np.linalg.norm(atom_a.coord - atom_b.coord)
            contacts.append({
                "chain_a": chain_a, "resseq_a": res_a.id[1],
                "resname_a": res_a.get_resname(), "atom_a": atom_a.get_name(),
                "chain_b": chain_b, "resseq_b": res_b.id[1],
                "resname_b": res_b.get_resname(), "atom_b": atom_b.get_name(),
                "distance": float(d),
            })

    # Build neighbor search on chain A for chain B interface
    ns_a = NeighborSearch(atoms_a)
    interface_b = set()
    for atom_b in atoms_b:
        neighbors = ns_a.search(atom_b.coord, cutoff, "A")
        for atom_a in neighbors:
            res_b = atom_b.get_parent()
            if res_b.id[0] != " ":
                continue
            interface_b.add((chain_b, res_b.id[1], res_b.get_resname()))

    # Sort
    interface_a = sorted(interface_a, key=lambda x: x[1])
    interface_b = sorted(interface_b, key=lambda x: x[1])

    result = {
        "chain_a": chain_a,
        "chain_b": chain_b,
        "cutoff": cutoff,
        "interface_a": interface_a,
        "interface_b": interface_b,
        "n_interface_a": len(interface_a),
        "n_interface_b": len(interface_b),
        "n_contacts": len(contacts),
    }

    print(f"Interface {chain_a}–{chain_b} (cutoff {cutoff}Å):")
    print(f"  Chain {chain_a}: {len(interface_a)} interface residues")
    print(f"  Chain {chain_b}: {len(interface_b)} interface residues")
    print(f"  Total atom-atom contacts: {len(contacts)}")

    return result, contacts
```

### 1.2 Classify Interface Contacts

```python
#!/opt/conda/envs/prot/bin/python
"""Classify interface contacts by type."""
import pandas as pd
import numpy as np

# Hydrogen bond donors and acceptors (simplified)
_HBOND_DONORS = {"N", "NE", "NH1", "NH2", "NZ", "ND1", "ND2", "NE1", "NE2", "OG", "OG1", "OH"}
_HBOND_ACCEPTORS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD"}
_CHARGED_POS = {"NZ", "NH1", "NH2", "NE"}  # Lys, Arg
_CHARGED_NEG = {"OD1", "OD2", "OE1", "OE2"}  # Asp, Glu
_AROMATIC_ATOMS = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH2", "CE3", "CZ2", "CZ3",
                   "NE1", "ND1", "CD2"}  # Phe, Tyr, Trp, His ring atoms
_HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}

def classify_contacts(contacts, hbond_cutoff=3.5, salt_cutoff=4.0,
                       aromatic_cutoff=5.5, hydrophobic_cutoff=4.5):
    """Classify each contact into type: hbond, salt_bridge, aromatic, hydrophobic, other.

    Args:
        contacts: list of contact dicts from extract_interface()
        *_cutoff: distance cutoffs for each contact type

    Returns:
        DataFrame with classified contacts + summary counts
    """
    classified = []
    for c in contacts:
        d = c["distance"]
        atom_a = c["atom_a"]
        atom_b = c["atom_b"]
        resname_a = c["resname_a"]
        resname_b = c["resname_b"]

        # Classify
        ctype = "other"

        # Salt bridge (charged-charged, < 4 Å)
        if d < salt_cutoff:
            if (atom_a in _CHARGED_POS and atom_b in _CHARGED_NEG) or \
               (atom_a in _CHARGED_NEG and atom_b in _CHARGED_POS):
                ctype = "salt_bridge"

        # H-bond (donor-acceptor, < 3.5 Å)
        if ctype == "other" and d < hbond_cutoff:
            if (atom_a in _HBOND_DONORS and atom_b in _HBOND_ACCEPTORS) or \
               (atom_a in _HBOND_ACCEPTORS and atom_b in _HBOND_DONORS):
                ctype = "hbond"

        # Aromatic (ring-ring, < 5.5 Å)
        if ctype == "other" and d < aromatic_cutoff:
            if atom_a in _AROMATIC_ATOMS and atom_b in _AROMATIC_ATOMS:
                ctype = "aromatic"

        # Hydrophobic (nonpolar-nonpolar, < 4.5 Å, carbon-carbon or hydrophobic res)
        if ctype == "other" and d < hydrophobic_cutoff:
            if (atom_a.startswith("C") and atom_b.startswith("C")) or \
               (resname_a in _HYDROPHOBIC_RES and resname_b in _HYDROPHOBIC_RES):
                ctype = "hydrophobic"

        row = dict(c)
        row["contact_type"] = ctype
        classified.append(row)

    df = pd.DataFrame(classified)

    # Summary
    if len(df) > 0:
        summary = df["contact_type"].value_counts().to_dict()
        print(f"Contact types: {summary}")
    else:
        summary = {}

    return df, summary
```

### 1.3 Interface Summary Report

```python
#!/opt/conda/envs/prot/bin/python
"""Generate interface summary report."""
import pandas as pd

def interface_summary(pdb_path, chain_a, chain_b, cutoff=5.0):
    """Full interface analysis: extraction + classification + per-residue contacts.

    Returns:
        summary dict, per-residue contact DataFrame, classified contacts DataFrame
    """
    iface, contacts = extract_interface(pdb_path, chain_a, chain_b, cutoff)
    classified_df, type_counts = classify_contacts(contacts)

    # Per-residue contact counts (on chain_b — the "design" chain)
    if len(classified_df) > 0:
        per_res = classified_df.groupby(["chain_b", "resseq_b", "resname_b"]).agg(
            n_contacts=("distance", "count"),
            n_hbond=("contact_type", lambda x: (x == "hbond").sum()),
            n_salt=("contact_type", lambda x: (x == "salt_bridge").sum()),
            n_aromatic=("contact_type", lambda x: (x == "aromatic").sum()),
            n_hydrophobic=("contact_type", lambda x: (x == "hydrophobic").sum()),
            min_distance=("distance", "min"),
        ).reset_index()
        per_res = per_res.sort_values("n_contacts", ascending=False)
    else:
        per_res = pd.DataFrame()

    summary = {
        **iface,
        "contact_types": type_counts,
        "has_salt_bridges": type_counts.get("salt_bridge", 0) > 0,
        "has_hbonds": type_counts.get("hbond", 0) > 0,
        "n_polar_contacts": type_counts.get("hbond", 0) + type_counts.get("salt_bridge", 0),
        "n_hydrophobic_contacts": type_counts.get("hydrophobic", 0),
    }

    return summary, per_res, classified_df
```

## Phase 2: Hotspot Identification

### 2.1 ESM-2 Alanine Scanning Proxy

```python
#!/opt/conda/envs/prot/bin/python
"""Identify interface hotspots using ESM-2 alanine scanning.

Logic: for each interface residue, compute the LLR of mutating it to Ala.
Positions where Ala is very unfavorable (large negative LLR from WT→Ala)
are positions where the wildtype residue is strongly preferred — likely
functionally/structurally important.

Conversely, positions where WT→Ala has near-zero LLR are tolerant of
mutation and safer to redesign.
"""
import torch
import pandas as pd
import numpy as np

def esm2_alanine_scan(sequence, interface_positions, model=None, alphabet=None,
                       batch_converter=None, device=None):
    """ESM-2 based alanine scanning of interface residues.

    For each interface position, computes:
    - WT log-likelihood at that position
    - Ala log-likelihood at that position
    - LLR (Ala - WT): negative = WT strongly preferred = potential hotspot

    Args:
        sequence: full amino acid sequence of the chain
        interface_positions: list of 0-indexed positions to scan
        model, alphabet, batch_converter, device: from load_esm2()

    Returns:
        DataFrame with per-position hotspot scores
    """
    # Load ESM-2 if not provided
    if model is None:
        import esm as esm_module
        model, alphabet = esm_module.pretrained.esm2_t33_650M_UR50D()
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

    ala_idx = alphabet.get_idx("A")

    rows = []
    for pos in interface_positions:
        wt_aa = sequence[pos]
        if wt_aa == "A":
            # Already Ala — no scan needed
            rows.append({
                "position": pos + 1,  # 1-indexed
                "wt_aa": wt_aa,
                "wt_ll": 0, "ala_ll": 0, "ala_llr": 0,
                "is_hotspot": False,
                "category": "ALREADY_ALA",
            })
            continue

        log_probs = torch.log_softmax(logits[pos + 1], dim=-1)  # +1 for BOS
        wt_idx = alphabet.get_idx(wt_aa)

        wt_ll = float(log_probs[wt_idx])
        ala_ll = float(log_probs[ala_idx])
        ala_llr = ala_ll - wt_ll  # negative means WT preferred over Ala

        # Classify
        # Large negative LLR: WT is strongly preferred → hotspot (mutation-sensitive)
        # Near-zero LLR: Ala is equally acceptable → tolerant position
        if ala_llr < -3.0:
            category = "STRONG_HOTSPOT"
            is_hotspot = True
        elif ala_llr < -1.5:
            category = "MODERATE_HOTSPOT"
            is_hotspot = True
        elif ala_llr < -0.5:
            category = "WEAK_HOTSPOT"
            is_hotspot = False
        else:
            category = "TOLERANT"
            is_hotspot = False

        rows.append({
            "position": pos + 1,
            "wt_aa": wt_aa,
            "wt_ll": float(wt_ll),
            "ala_ll": float(ala_ll),
            "ala_llr": float(ala_llr),
            "is_hotspot": is_hotspot,
            "category": category,
        })

    df = pd.DataFrame(rows)

    if len(df) > 0:
        n_hot = df["is_hotspot"].sum()
        print(f"Alanine scan: {n_hot}/{len(df)} interface positions are hotspots")
        print(f"  STRONG: {(df['category']=='STRONG_HOTSPOT').sum()}, "
              f"MODERATE: {(df['category']=='MODERATE_HOTSPOT').sum()}, "
              f"WEAK: {(df['category']=='WEAK_HOTSPOT').sum()}, "
              f"TOLERANT: {(df['category']=='TOLERANT').sum()}")

    return df
```

### 2.2 Hotspot Classification Thresholds

| Ala LLR (WT→Ala) | Category | Interpretation | Design action |
|-------------------|----------|---------------|---------------|
| < -3.0 | STRONG_HOTSPOT | WT residue strongly required | **FREEZE** — do not mutate |
| -3.0 to -1.5 | MODERATE_HOTSPOT | WT residue preferred | **CAUTION** — mutate only with clear rationale |
| -1.5 to -0.5 | WEAK_HOTSPOT | Some preference for WT | OK to explore conservative mutations |
| > -0.5 | TOLERANT | Ala equally acceptable | **SAFE to redesign** |

### 2.3 Combining Hotspots with Contact Type

```python
#!/opt/conda/envs/prot/bin/python
"""Combine hotspot scores with contact classification."""
import pandas as pd

def hotspot_contact_profile(hotspot_df, per_residue_contacts_df):
    """Merge hotspot scores with per-residue contact information.

    This gives the complete picture: which positions are energetically
    important AND what types of contacts they make.

    Args:
        hotspot_df: from esm2_alanine_scan()
        per_residue_contacts_df: from interface_summary()

    Returns:
        merged DataFrame with hotspot + contact info per residue
    """
    # Merge on position
    merged = hotspot_df.merge(
        per_residue_contacts_df,
        left_on="position",
        right_on="resseq_b",
        how="left",
    )

    # Fill missing contact counts with 0
    for col in ["n_contacts", "n_hbond", "n_salt", "n_aromatic", "n_hydrophobic"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    # Design recommendation
    recommendations = []
    for _, row in merged.iterrows():
        if row["category"] == "STRONG_HOTSPOT":
            rec = "FREEZE"
        elif row["category"] == "MODERATE_HOTSPOT":
            if row.get("n_hbond", 0) > 0 or row.get("n_salt", 0) > 0:
                rec = "FREEZE — polar contacts"
            else:
                rec = "CAUTION — conservative mutations only"
        elif row["category"] == "WEAK_HOTSPOT":
            rec = "EXPLORE — conservative mutations"
        elif row["category"] == "ALREADY_ALA":
            rec = "EXPLORE — already Ala"
        else:
            rec = "REDESIGN — tolerant position"
        recommendations.append(rec)

    merged["design_recommendation"] = recommendations

    return merged
```

### 2.4 Burial Analysis (SASA-Based Core vs Surface vs Interface Classification)

```python
#!/opt/conda/envs/prot/bin/python
"""Classify residues as CORE / SURFACE / INTERFACE based on SASA.

This is critical for interpreting ESM-2 LLR results:
- CORE mutations with low LLR → likely truly destabilizing (trust ESM-2)
- INTERFACE mutations with low LLR → may be ESM-2 evolutionary bias (cross-check)
- SURFACE mutations with low LLR → usually safe to discount
"""
from Bio.PDB import PDBParser, DSSP, ShrakeRupley
import numpy as np
import pandas as pd

def compute_residue_sasa(pdb_path, chain_id=None):
    """Compute per-residue SASA using Biopython ShrakeRupley.

    Args:
        pdb_path: PDB file
        chain_id: specific chain (None = all chains)

    Returns:
        DataFrame with resseq, resname, sasa, and burial classification
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    model = list(structure.get_models())[0]

    # Compute SASA
    sr = ShrakeRupley()
    sr.compute(model, level="R")  # residue-level SASA

    rows = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            sasa = residue.sasa
            rows.append({
                "chain": chain.id,
                "resseq": residue.id[1],
                "resname": residue.get_resname(),
                "sasa": float(sasa),
            })

    df = pd.DataFrame(rows)

    if len(df) > 0:
        # Classify burial
        # Thresholds based on standard residue max SASA (~200-250 Å² for large residues)
        # Relative SASA would be better, but absolute works for classification
        df["burial"] = df["sasa"].apply(_classify_burial)

    return df

def _classify_burial(sasa):
    """Classify residue burial from absolute SASA.

    Thresholds (approximate, for standard amino acids):
    - CORE: SASA < 10 Å² (buried, minimal solvent exposure)
    - PARTIAL: 10-40 Å² (partially buried)
    - SURFACE: > 40 Å² (exposed)
    """
    if sasa < 10:
        return "CORE"
    elif sasa <= 40:
        return "PARTIAL"
    else:
        return "SURFACE"

def compute_relative_sasa(pdb_path, chain_id=None):
    """Compute relative SASA (fraction of max for each residue type).

    More accurate than absolute SASA for burial classification.
    Uses Gly-X-Gly reference values.
    """
    # Max SASA reference values (Tien et al. 2013, Gly-X-Gly tripeptide)
    MAX_SASA = {
        "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
        "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
        "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
        "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
    }

    df = compute_residue_sasa(pdb_path, chain_id)

    if len(df) > 0:
        df["max_sasa"] = df["resname"].map(MAX_SASA).fillna(200)
        df["relative_sasa"] = df["sasa"] / df["max_sasa"]
        df["burial_rel"] = df["relative_sasa"].apply(_classify_burial_relative)

    return df

def _classify_burial_relative(rel_sasa):
    """Classify burial from relative SASA (more accurate).

    - CORE: relSASA < 0.05 (< 5% exposed)
    - PARTIAL: 0.05-0.20
    - SURFACE: > 0.20 (> 20% exposed)
    """
    if rel_sasa < 0.05:
        return "CORE"
    elif rel_sasa <= 0.20:
        return "PARTIAL"
    else:
        return "SURFACE"


def annotate_interface_burial(pdb_path, interface_residues, chain_id):
    """Add burial classification to interface residues.

    This is the key diagnostic: an interface residue that is CORE-like
    (low SASA even in the complex) is structurally critical and ESM-2
    LLR rejection should be trusted. An interface residue that is
    SURFACE/PARTIAL in the unbound state is more likely to tolerate
    mutation — ESM-2 LLR may be overly conservative.

    Args:
        pdb_path: PDB of the complex
        interface_residues: list of (chain, resseq, resname) tuples
        chain_id: chain to analyze

    Returns:
        DataFrame with burial annotation per interface residue
    """
    sasa_df = compute_relative_sasa(pdb_path, chain_id)

    # Map interface residues
    iface_set = set((r[1] if len(r) > 1 else r) for r in interface_residues)

    # Filter to interface
    iface_sasa = sasa_df[sasa_df["resseq"].isin(iface_set)].copy()

    if len(iface_sasa) > 0:
        n_core = (iface_sasa["burial_rel"] == "CORE").sum()
        n_partial = (iface_sasa["burial_rel"] == "PARTIAL").sum()
        n_surface = (iface_sasa["burial_rel"] == "SURFACE").sum()
        print(f"Interface burial (chain {chain_id}): "
              f"CORE={n_core}, PARTIAL={n_partial}, SURFACE={n_surface}")

    return iface_sasa
```

### 2.5 Interpreting ESM-2 LLR with Burial Context

| Burial | ESM-2 LLR < -2 | Action |
|--------|----------------|--------|
| **CORE** | Trust ESM-2 — mutation likely destabilizing | FREEZE or use physics (FoldX) to confirm |
| **PARTIAL** | Ambiguous — could be fold or function | Flag for review; consider conservative mutations |
| **SURFACE** | Likely ESM-2 evolutionary bias | Discount LLR; trust fold-back + interface contacts instead |

> **This is why 1BRS A36S was rejected by LLR but passed fold-back + interface gates:** if A36 is CORE/PARTIAL in the complex, ESM-2 is probably right. If it's SURFACE at the interface, ESM-2 is likely over-penalizing a legitimate interface mutation.

## Phase 3: Design Constraint Generation

### 3.1 Generate prot-seqdesign Constraints from Interface Analysis

```python
#!/opt/conda/envs/prot/bin/python
"""Convert interface analysis into prot-seqdesign constraints."""

def generate_design_constraints(hotspot_profile_df, chain_id,
                                 freeze_strong=True, freeze_moderate_polar=True):
    """Generate residue lists for prot-seqdesign from hotspot analysis.

    Outputs:
    - fixed_residues: positions to FREEZE (don't redesign)
    - redesign_residues: positions safe to redesign
    - caution_residues: positions that need careful consideration

    Args:
        hotspot_profile_df: from hotspot_contact_profile()
        chain_id: chain ID for formatting
        freeze_strong: freeze STRONG_HOTSPOT positions (default True)
        freeze_moderate_polar: freeze MODERATE_HOTSPOT with polar contacts (default True)

    Returns:
        dict with residue lists formatted for LigandMPNN --redesigned_residues
    """
    fixed = []
    redesign = []
    caution = []

    for _, row in hotspot_profile_df.iterrows():
        pos = int(row["position"])
        label = f"{chain_id}{pos}"
        cat = row["category"]

        if cat == "STRONG_HOTSPOT" and freeze_strong:
            fixed.append(label)
        elif cat == "MODERATE_HOTSPOT":
            if freeze_moderate_polar and row.get("n_hbond", 0) + row.get("n_salt", 0) > 0:
                fixed.append(label)
            else:
                caution.append(label)
        elif cat == "ALREADY_ALA":
            redesign.append(label)
        elif cat in ("WEAK_HOTSPOT", "TOLERANT"):
            redesign.append(label)
        else:
            redesign.append(label)

    result = {
        "fixed_residues": fixed,
        "redesign_residues": redesign,
        "caution_residues": caution,
        "redesign_string": " ".join(redesign),  # for --redesigned_residues
        "fixed_string": " ".join(fixed),
    }

    print(f"Design constraints for chain {chain_id}:")
    print(f"  FREEZE ({len(fixed)}): {' '.join(fixed)}")
    print(f"  REDESIGN ({len(redesign)}): {' '.join(redesign)}")
    print(f"  CAUTION ({len(caution)}): {' '.join(caution)}")
    print(f"\nFor prot-seqdesign --redesigned_residues: \"{result['redesign_string']}\"")

    return result
```

### 3.2 Decision: Full Interface Redesign vs Targeted Mutations

```
You want to optimize an interface. Which approach?
│
├── Goal: improve affinity for SAME target
│   ├── Hotspot analysis shows >50% tolerant positions?
│   │   └── YES → Full interface redesign (redesign all tolerant + caution positions)
│   │   └── NO → Targeted mutations only (pick 2-3 tolerant positions, leave rest fixed)
│   └── Use prot-seqdesign with --redesigned_residues from Phase 3.1
│
├── Goal: change specificity (bind different target)
│   ├── This is harder — need to DISRUPT old contacts AND BUILD new ones
│   ├── Identify hotspots for OLD target → these must change
│   ├── Model new target interface → identify positions that need new contacts
│   └── Usually requires backbone redesign (prot-backbone-gen, future skill)
│
├── Goal: disrupt binding (make non-binder)
│   ├── Mutate STRONG_HOTSPOT positions to Ala or opposite charge
│   ├── 2-3 hotspot mutations usually sufficient
│   └── Verify with prot-stability that fold is maintained
│
└── Goal: improve stability without losing binding
    ├── ONLY redesign non-interface positions (surface away from interface)
    ├── Use prot-seqdesign with interface positions FROZEN
    └── Verify interface contacts preserved in designed structure (prot-structure fold-back)
```

## Phase 4: Interface Quality Assessment for Designs

### 4.1 Compare Designed vs Wildtype Interface

```python
#!/opt/conda/envs/prot/bin/python
"""Compare interface properties between wildtype and designed complex."""

def compare_interfaces(wt_pdb, design_pdb, chain_a, chain_b, cutoff=5.0):
    """Compare interface contacts between wildtype and designed structures.

    Use after fold-back validation (prot-structure) to check if
    designed mutations preserve or improve interface contacts.

    Args:
        wt_pdb: wildtype complex PDB
        design_pdb: designed complex PDB (from fold-back or modeling)
        chain_a, chain_b: chain IDs
        cutoff: contact distance cutoff

    Returns:
        comparison dict with gained/lost/maintained contacts
    """
    wt_summary, wt_per_res, wt_contacts = interface_summary(wt_pdb, chain_a, chain_b, cutoff)
    des_summary, des_per_res, des_contacts = interface_summary(design_pdb, chain_a, chain_b, cutoff)

    # Compare contact counts by type
    comparison = {
        "wt_n_contacts": wt_summary["n_contacts"],
        "des_n_contacts": des_summary["n_contacts"],
        "wt_n_polar": wt_summary["n_polar_contacts"],
        "des_n_polar": des_summary["n_polar_contacts"],
        "wt_n_hydrophobic": wt_summary["n_hydrophobic_contacts"],
        "des_n_hydrophobic": des_summary["n_hydrophobic_contacts"],
    }

    # Compare interface residue sets
    wt_res_b = set((r[1], r[2]) for r in wt_summary["interface_b"])
    des_res_b = set((r[1], r[2]) for r in des_summary["interface_b"])

    comparison["maintained_residues"] = len(wt_res_b & des_res_b)
    comparison["lost_residues"] = len(wt_res_b - des_res_b)
    comparison["gained_residues"] = len(des_res_b - wt_res_b)

    # Assessment
    lost_polar = wt_summary["n_polar_contacts"] - des_summary.get("n_polar_contacts", 0)
    if comparison["lost_residues"] == 0 and lost_polar <= 0:
        comparison["assessment"] = "PRESERVED"
        comparison["action"] = "Interface contacts maintained or improved."
    elif comparison["lost_residues"] <= 2 and lost_polar <= 1:
        comparison["assessment"] = "MINOR_CHANGE"
        comparison["action"] = "Small interface changes. Check if lost contacts are critical."
    else:
        comparison["assessment"] = "SIGNIFICANT_CHANGE"
        comparison["action"] = (f"Lost {comparison['lost_residues']} interface residues "
                               f"and {lost_polar} polar contacts. Review carefully.")

    print(f"\nInterface comparison (WT vs Design):")
    print(f"  Contacts: {comparison['wt_n_contacts']} → {comparison['des_n_contacts']}")
    print(f"  Polar: {comparison['wt_n_polar']} → {comparison['des_n_polar']}")
    print(f"  Interface residues: {comparison['maintained_residues']} maintained, "
          f"{comparison['lost_residues']} lost, {comparison['gained_residues']} gained")
    print(f"  Assessment: {comparison['assessment']}")

    return comparison
```

## Phase 5: Interface Decision Tree — Full Workflow

```
You have a protein complex and need to analyze/design at the interface.
│
├── Step 1: Get the structure
│   ├── Experimental complex PDB? → Use it (best)
│   ├── Only individual structures? → Need docking or AF2-multimer (not covered here)
│   └── Gate: structure quality check (prot-structure Phase 2) BEFORE interface analysis
│
├── Step 2: Map the interface (Phase 1)
│   ├── extract_interface() with cutoff=5Å
│   ├── classify_contacts() → get contact types
│   └── interface_summary() → per-residue contact profile
│
├── Step 3: Identify hotspots (Phase 2)
│   ├── esm2_alanine_scan() on interface residues
│   ├── hotspot_contact_profile() → merge with contacts
│   └── Review: STRONG/MODERATE hotspots typically 2-5 residues
│
├── Step 4: Generate design constraints (Phase 3)
│   ├── generate_design_constraints() → fixed/redesign/caution lists
│   └── Feed --redesigned_residues into prot-seqdesign
│
├── Step 5: Design sequences (prot-seqdesign)
│   └── Use interface constraints from Step 4
│
├── Step 6: Score designs
│   ├── prot-stability: ESM-2 LLR for fold stability
│   ├── prot-structure: fold-back validation
│   └── prot-interface: compare_interfaces() for binding preservation
│
└── Step 7: Final candidate selection
    └── Must pass ALL three: stability OK + fold-back OK + interface preserved
```

## Phase 6: Hard Gates for Interface Design

```python
#!/opt/conda/envs/prot/bin/python
"""Interface-specific hard gates for design screening."""

def interface_hard_gates(hotspot_profile, design_mutations):
    """Check if design mutations violate interface constraints.

    Args:
        hotspot_profile: from hotspot_contact_profile()
        design_mutations: list of (position_1indexed, wt_aa, mut_aa)

    Returns:
        pass/fail, list of violations
    """
    violations = []
    warnings = []

    hotspot_dict = {}
    for _, row in hotspot_profile.iterrows():
        hotspot_dict[int(row["position"])] = row

    for pos, wt, mut in design_mutations:
        if pos not in hotspot_dict:
            continue  # not an interface position

        row = hotspot_dict[pos]
        cat = row["category"]

        # Gate 1: Never mutate STRONG_HOTSPOT
        if cat == "STRONG_HOTSPOT":
            violations.append(
                f"VIOLATION: {wt}{pos}{mut} mutates STRONG_HOTSPOT "
                f"(ala_LLR={row['ala_llr']:.2f})")

        # Gate 2: MODERATE_HOTSPOT with polar contacts — warn
        elif cat == "MODERATE_HOTSPOT":
            has_polar = row.get("n_hbond", 0) + row.get("n_salt", 0) > 0
            if has_polar:
                violations.append(
                    f"VIOLATION: {wt}{pos}{mut} mutates MODERATE_HOTSPOT "
                    f"with polar contacts")
            else:
                warnings.append(
                    f"WARNING: {wt}{pos}{mut} mutates MODERATE_HOTSPOT "
                    f"(no polar contacts — may be OK)")

    passed = len(violations) == 0

    return {
        "passed": passed,
        "n_violations": len(violations),
        "n_warnings": len(warnings),
        "violations": violations,
        "warnings": warnings,
    }
```

## Phase 7: VRAM Budget & Timing

| Task | VRAM | Wall time | Notes |
|------|------|-----------|-------|
| Interface extraction | 0 GB | ~1-2 sec | CPU, Biopython |
| Contact classification | 0 GB | ~1 sec | CPU, pure Python |
| ESM-2 alanine scan (load) | ~2 GB | ~10 sec | One-time model load |
| ESM-2 alanine scan (scan) | ~2 GB | ~1 sec | Single forward pass |
| Interface comparison | 0 GB | ~2-3 sec | CPU, two PDBs |

> ESM-2 can stay loaded alongside LigandMPNN on 12 GB VRAM.

## Failure Modes

1. **Designing without interface analysis.** The #1 failure mode. Blind redesign of all interface residues will likely disrupt binding. Always map hotspots first and freeze critical positions.

2. **Confusing fold stability with binding importance.** ESM-2 alanine scanning captures evolutionary importance, which includes both fold stability and function. A STRONG_HOTSPOT might be fold-critical (buried Trp in the core) rather than binding-critical (interfacial contact). Cross-reference with contact type to distinguish: if it's a STRONG_HOTSPOT WITH polar interface contacts, it's likely binding-critical.

3. **Trusting a single cutoff.** A 5Å cutoff gives a different interface than 4Å or 8Å. For design constraints, 5Å is standard. For "core interface" contacts, tighten to 4Å. For "extended interface" effects, expand to 8Å.

4. **Ignoring water-mediated contacts.** The current implementation only detects direct contacts. Water-mediated hydrogen bonds can be critical but are invisible here. If a key H-bond seems missing from the contact map, check for bridging waters in the crystal structure.

5. **Over-constraining the design.** Freezing too many positions leaves too few to redesign, reducing design diversity. If >70% of interface positions are hotspots, the interface is already highly optimized — consider designing away from the interface instead.

6. **Short peptide interfaces.** For peptide-protein interfaces (<30 aa peptide), the peptide may have very few interface residues (5-10), making each position more precious. Be more conservative with redesign — even "TOLERANT" positions matter when there are only a few.

## One-Sentence Rule

**Map the interface and identify hotspots before designing — freeze strong hotspots, redesign tolerant positions, and always verify that designed mutations preserve critical contacts.**
