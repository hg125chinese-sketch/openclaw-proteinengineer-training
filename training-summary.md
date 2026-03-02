# ProteinEngineer Training Summary

Date: 2026-03-01
Trainer: Claude (Anthropic)
Agent: ProteinEngineer (PR)
Platform: OpenClaw

---

## 1. Training Overview

### Methodology
7-step training framework adapted from ChemicalExpert:
1. Pre-assessment (14 questions, 4 dimensions)
2. Deep Probing (5 targeted questions → DP score)
3. Gap Analysis (tool vs methodology gaps)
4. Environment Setup (conda, tools, dependencies)
5. Skill Design + Guided Practice
6. Behavioral Correction (skill iteration based on real-run feedback)
7. Validation (multi-skill integration tests)

### Pre-Assessment Results
- Methodology: 9/10 (strong: knows design principles, DBTL logic)
- Tool Competence: 4/10 (weak: never ran LigandMPNN, ESMFold, ESM-2 in this environment)
- Decision Making: 7/10 (good judgment, needs calibrated thresholds)
- Integration: 5/10 (understands pipeline concept, lacks operational workflow)

### Deep Probe Score
- DP2 (after probe): 2/10 — critical tool gap confirmed

---

## 2. Skills Developed (7 total)

| # | Skill | Version | Description | Key Iteration |
|---|-------|---------|-------------|---------------|
| 1 | prot-seqdesign | v2 | LigandMPNN + ESM-IF sequence design | v1→v2: CLI auto-detect, ESM-IF signature fix, biotite patch, dual output parser |
| 2 | prot-structure | **v4** | ESMFold structure prediction & fold-back | v1→v2: API-first (openfold blocker); v3: pLDDT unit fix, short peptide RMSD, bonded exclusion; **v4: WT sanity check** |
| 3 | prot-stability | v2 | ESM-2 LLR stability scoring | v1→v2: interface-mutation mode with burial-aware gates |
| 4 | prot-interface | v2 | Interface analysis, hotspots, design constraints | v1→v2: SASA burial analysis (CORE/PARTIAL/SURFACE) |
| 5 | prot-dbtl | v1 | DBTL cycle orchestration | — |
| 6 | prot-msa | v1 | Conservation analysis + LLR calibration | — |
| 7 | prot-developability | v1 | Red-light screening before experiments | — |

### One-Sentence Rules

1. **prot-seqdesign**: Design on validated structures, compare every score to wildtype, check your score direction (NEW vs LEGACY), and cross-validate LigandMPNN with ESM-IF before advancing.
2. **prot-structure**: Gate structure quality before design, run WT sanity check before fold-back, and if ESMFold can't fold the wildtype — downgrade to diagnostic.
3. **prot-stability**: Screen with ESM-2 LLR, use interface mode with burial context for interface tasks, compare to wildtype, diagnose disagreements rather than averaging.
4. **prot-interface**: Map interface and hotspots before designing — freeze strong hotspots, redesign tolerant positions, verify contacts preserved after design.
5. **prot-dbtl**: Plan every cycle before executing, log every decision with rationale, let experimental data override computational predictions.
6. **prot-msa**: Use family conservation to calibrate LLR — variable positions tolerate mutations that LLR penalizes, conserved positions are dangerous even when LLR is mild.
7. **prot-developability**: Run every candidate through the red-light table before ordering DNA — free cysteines, deamidation, and aggregation kill more designs than bad scores.

---

## 3. Practice Runs

### 3.1 1YCR (p53-MDM2) — Skills 1-4 individual validation

**Target**: 13aa p53 peptide (chain B) bound to MDM2 (chain A)

**Key results**:
- prot-seqdesign: 50 designs → 4 unique, dominant ETFEDLWKKLPQS (46/50)
- prot-structure v3: fold-back PASS (RMSD<2Å, short peptide mode)
- prot-stability: all Top 3 REJECT (mean_LLR -3.28 to -3.51)
- prot-interface: 10/11 interface positions STRONG_HOTSPOT → entire peptide frozen

**Key learnings**:
- Short peptide edge cases exposed in all 4 skills
- ESM-2 LLR magnitudes inflated for short sequences
- TM-score unreliable for <30aa (switched to RMSD)
- API B-factors in 0-1 range (auto-detect + ×100 fix)
- Clash detection needed bonded exclusion (min_seq_sep=3)

### 3.2 1BRS (barnase-barstar) — 4-skill integration test

**Target**: barstar chain D (89aa) interface with barnase chain A (110aa)

**Key results**:
- Interface: 17/19 STRONG_HOTSPOT (89.5%) → only 2 redesign positions
- seqdesign: 1 unique candidate (A36S) — too few free positions
- stability: mean_LLR=-6.19 → catastrophic FAIL (standard mode)
- structure: fold-back TM=0.97, pLDDT=87 → PASS
- interface: mutation in redesign set → PASS

**Critical finding**: fold-back PASS + interface PASS + stability FAIL
→ Diagnosed as ESM-2 evolutionary bias against interface mutations
→ Led to prot-stability v2 (interface-mutation mode) and prot-interface v2 (burial analysis)

### 3.3 1Z92 (IL-2 / IL-2Rα) — 7-skill end-to-end validation

**Target**: IL-2 chain A (121aa) interface with IL-2Rα chain B (123aa)

**Key results**:
- Interface: 15/22 STRONG_HOTSPOT (68.2%)
- MSA: 0 consensus mutations suggested (WT already at ESM-2 consensus)
- Constraints: tight-mode → 2 redesign positions
- seqdesign: 1 unique candidate (F42Y)
- stability (interface mode): mean_LLR=-3.33 → PASS (relaxed gate)
- structure: WT pLDDT≈45, TM≈0.49 → ESMFold unreliable for IL-2
- developability: ODD_CYSTEINES(3), instability>40, 1 deamidation — same as WT

**Critical finding**: ESMFold fold-back gate falsely rejects everything because WT itself scores poorly (disulfide-rich protein)
→ Led to prot-structure v4 (WT sanity check: HARD_GATE / DIAGNOSTIC / SKIP)

**Secondary finding**: developability needs WT-relative mode (don't penalize designs for issues the WT already has)

---

## 4. Skill Iteration History

```
prot-seqdesign:
  v1 (initial) → v2 (5 fixes from 1YCR practice)

prot-structure:
  v1 (local ESMFold) → v2 (API-first, openfold blocker)
  → v3 (pLDDT unit fix, short peptide RMSD, clash bonded exclusion)
  → v4 (WT sanity check from 1Z92 debug)

prot-stability:
  v1 (standard gates) → v2 (interface-mutation mode + burial gates from 1BRS)

prot-interface:
  v1 (hotspot + constraints) → v2 (SASA burial analysis from 1BRS diagnosis)

prot-dbtl: v1 (no iteration needed)
prot-msa: v1 (no iteration needed)
prot-developability: v1 (no iteration needed)
```

Every iteration was driven by **real practice run failures**, not hypothetical improvements.

---

## 5. Pipeline Architecture

```
prot-interface (Step 2)          prot-msa (Step 3)
  ├─ hotspot analysis              ├─ conservation
  ├─ burial/SASA                   └─ LLR calibration
  └─ design constraints                    │
         │                                 │
         ▼                                 ▼
prot-seqdesign (Step 5)     prot-stability (Step 6)
  ├─ LigandMPNN                 ├─ ESM-2 LLR
  └─ ESM-IF cross-val           ├─ interface mode + burial
                                └─ conservation calibration
         │                                 │
         ▼                                 ▼
prot-structure (Step 7)     prot-developability (Step 9)
  ├─ WT sanity check             ├─ chemical liabilities
  ├─ fold-back validation        ├─ aggregation/solubility
  └─ HARD/DIAGNOSTIC/SKIP       └─ expression compatibility
         │                                 │
         └──────────────┬──────────────────┘
                        ▼
               prot-dbtl (orchestration)
                 ├─ cycle planning
                 ├─ composite scoring
                 ├─ gate summary
                 ├─ decision logging
                 └─ experimental feedback → next cycle
```

---

## 6. Key Methodological Insights

### What worked in the training methodology
1. **Practice-first skill iteration**: every skill version was driven by real run failures, not speculation
2. **Incremental complexity**: 1YCR (short peptide, single skill) → 1BRS (4-skill) → 1Z92 (7-skill)
3. **PR's autonomous diagnostics**: agent independently identified root causes (openfold blocker, pLDDT units, TM-score unreliability) and wrote structured feedback
4. **skill-feedback.md as learning journal**: captures what worked, what broke, and why — reusable across training sessions

### Recurring patterns across practice runs
1. **Standard thresholds break on edge cases**: short peptides, disulfide-rich proteins, tight interfaces — every practice run exposed a new edge case that required adaptive gate logic
2. **ESM-2 is a global prior, not ground truth**: penalizes legitimate interface mutations, inflates LLR for short sequences, gives false confidence for conserved-but-wrong mutations
3. **WT baseline comparison is essential**: for stability (ΔLLR not absolute), developability (same flags as WT = OK), and fold-back (if WT fails, designs will too)

### Known remaining gaps
1. **No physics-based ΔΔG** (FoldX/Rosetta require academic licenses)
2. **ESMFold weak on disulfide-rich / multi-domain proteins** (mitigated by WT sanity check, not solved)
3. **Tight interfaces (>50% hotspot) produce too few redesign positions** (need backbone redesign or relaxed constraints)
4. **Developability needs WT-relative mode** (flagged but not yet implemented)
5. **No complex structure prediction** (single-chain fold-back only, no AF2-multimer)

---

## 7. Files & Artifacts

### Workspace
- PLAYBOOK.md: operational rules for all 7 skills
- skill-review.md: PR's self-assessment (updated per practice run)
- skill-feedback.md: detailed per-run feedback with fix suggestions

### Skills (vault + workspace)
```
skills/heng/
├── prot-seqdesign/SKILL.md    (v2)
├── prot-structure/SKILL.md    (v4)
├── prot-stability/SKILL.md    (v2)
├── prot-interface/SKILL.md    (v2)
├── prot-dbtl/SKILL.md         (v1)
├── prot-msa/SKILL.md          (v1)
└── prot-developability/SKILL.md (v1)
```

### Practice run directories
```
runs/
├── 1YCR_seqdesign/        (prot-seqdesign v1→v2)
├── 1YCR_structure_v3/     (prot-structure v2→v3 validation)
├── 1YCR_stability/        (prot-stability v1)
├── 1YCR_interface/        (prot-interface v1)
├── 1BRS_integration/      (4-skill integration → stability v2, interface v2)
└── 1Z92_cycle1/           (7-skill end-to-end → structure v4)
    └── debug_foldback/    (WT sanity check development)
```

### GitHub
- Repository: hg125chinese-sketch/openclaw-chemicalexpert-training
- Training guide release: pending
