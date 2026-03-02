# Teaching AI Agents to Engineer Proteins: A Systematic Training Methodology

How I trained an LLM agent to autonomously run protein design workflows — and what broke along the way taught us more than what worked.

## TL;DR

I trained an AI agent ("ProteinEngineer") with 7 specialized protein engineering skills using a systematic 7-step methodology adapted from an earlier drug discovery agent. The agent learned to autonomously plan and execute interface redesign campaigns: analyzing binding interfaces, generating constrained designs, scoring stability, validating folds, and screening for manufacturability.

When tested on three progressively harder targets (13aa peptide → 89aa protein → 121aa cytokine), each practice run exposed edge cases that drove skill iterations — short peptide gate failures, ESM-2 evolutionary bias against interface mutations, and ESMFold's inability to fold disulfide-rich proteins.

This post describes the methodology, the skills, and what happened when the agent's tools disagreed with each other.

## The Problem

Protein language models know a lot about sequences. Structure predictors can fold most proteins. But orchestrating a complete design campaign — from interface analysis through sequence design to experimental readiness — requires judgment calls that no single tool provides.

I wanted an agent that could:

- Run protein interface redesign end-to-end (structure → interface analysis → constrained design → multi-tool evaluation)
- Know when its tools are reliable and when they aren't
- Diagnose conflicts between tools rather than blindly trusting any one of them
- Adapt its gate thresholds to the specific protein being designed

The platform is OpenClaw, an open-source agent framework. The underlying model is Claude via Anthropic's API. The agent runs in a Docker container with GPU support on NixOS/WSL2, with access to conda environments, LigandMPNN, ESMFold, ESM-2, and fair-esm.

## The 7-Step Training Methodology

Adapted from ChemicalExpert training (see [openclaw-chemicalexpert-training](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training)), with protein-specific adjustments:

**Step 1: Pre-Assessment.** 14 diagnostic questions across four dimensions: methodology, tool competence, decision making, and integration capability. PE scored 9/10 on methodology but 4/10 on tool competence — knew the theory but had never run the tools.

**Step 2: Deep Probing.** 5 targeted questions exposing the gap between "knows about LigandMPNN" and "can run LigandMPNN in this container with this fork's CLI flags." Deep probe score: 2/10.

**Step 3: Gap Analysis.** Tool gaps dominated. The agent understood protein design principles but couldn't execute any of them. Training focused entirely on operational capability.

**Step 4: Skill Design.** Each skill is a structured document (200-600 lines) with Phase 0 environment verification, executable code blocks, decision trees, hard gates, failure modes, and a one-sentence rule.

**Step 5: Guided Practice.** Real targets, not toy problems. Every practice run used actual PDB structures with known biology, producing results that could be evaluated against experimental literature.

**Step 6: Behavioral Correction.** When skills produced wrong results, the agent was challenged to diagnose the root cause. Every skill iteration was driven by a real failure, not hypothetical improvements.

**Step 7: Validation.** Multi-skill integration tests on progressively harder targets, culminating in a 7-skill end-to-end run.

## The 7 Skills

| # | Skill | Version | Core Capability |
|---|-------|---------|----------------|
| 1 | prot-seqdesign | v2 | LigandMPNN sequence design + ESM-IF cross-validation |
| 2 | prot-structure | v4 | ESMFold structure prediction, fold-back validation, WT sanity check |
| 3 | prot-stability | v2 | ESM-2 LLR scoring with interface-mutation mode + burial-aware gates |
| 4 | prot-interface | v2 | Interface analysis, hotspot identification, SASA burial classification |
| 5 | prot-msa | v1 | Conservation analysis, consensus design, LLR calibration |
| 6 | prot-developability | v2 | Chemical liabilities, aggregation, expression compatibility, WT-relative mode |
| 7 | prot-dbtl | v1 | Design-Build-Test-Learn cycle orchestration |

Each skill document is in `skills/*/SKILL.md`.

## Skill Design Principles

**Philosophy over procedure.** "Use family conservation to calibrate ESM-2 LLR" beats 50 lines of API documentation.

**Decision trees over checklists.** "If WT fold-back pLDDT < 60, downgrade to diagnostic mode" matches real engineering decisions.

**Failure modes are non-negotiable.** Every skill includes what failure looks like and what causes it. These sections were written from real failures, not imagination.

**Gates, not guidelines.** Hard thresholds (TM ≥ 0.85 for PASS, mean_LLR < -2.0 for REJECT) that must be passed. But gates adapt: interface mode relaxes LLR gates when fold-back and interface gates are also applied.

**One sentence per skill.** If you can't state the core rule in one sentence, the skill is too broad.

**WT baseline comparison.** Every metric should be compared to the wildtype. An LLR of -3 means nothing without knowing what the WT scores.

## The Practice Runs

### Run 1: 1YCR (p53-MDM2) — Short Peptide Edge Cases

A 13-amino-acid p53 peptide bound to MDM2. Chosen as a simple first target. Exposed edge cases in three skills simultaneously.

**What broke:**

- ESMFold API returns B-factors in 0-1 range, not 0-100 → pLDDT gate falsely rejected everything
- Clash detection counted bonded neighbors as clashes → 11 false clashes on a 13aa peptide
- TM-score normalization is statistically unreliable for sequences <30aa → valid designs scored TM=0.35

**What was fixed:**

- prot-structure v3: auto-detect pLDDT units (if max ≤ 1.0, multiply by 100)
- prot-structure v3: bonded exclusion with min_seq_sep=3
- prot-structure v3: short peptides use RMSD as primary metric, not TM-score

### Run 2: 1BRS (Barnase-Barstar) — Tool Disagreement

Standard-sized proteins (110aa + 89aa). First 4-skill integration test. Exposed a fundamental conflict between tools.

**What happened:**

- Interface analysis found 89.5% STRONG hotspots → only 2 positions available for redesign
- LigandMPNN produced 1 unique design (A36S) — almost no sequence diversity
- Fold-back: TM = 0.97, pLDDT = 87 → **PASS**
- Interface gate: mutation in allowed set → **PASS**
- Stability: mean_LLR = -6.19 → **catastrophic FAIL**

**The diagnosis:**

ESM-2 LLR captures evolutionary conservation, not just thermodynamic stability. Interface positions are evolutionarily conserved because they're functionally important — but that's exactly where you want to mutate for affinity maturation.

The agent correctly identified this as "ESM-2 evolutionary bias" rather than a true stability problem.

**What was built:**

- prot-stability v2: interface-mutation mode with relaxed gates (mean_LLR -2→-4, catastrophic -4→-6)
- prot-interface v2: SASA-based burial analysis to distinguish CORE (trust ESM-2) from SURFACE (discount ESM-2)
- Integration: burial context feeds into stability gates — CORE mutations still use strict thresholds

### Run 3: 1Z92 (IL-2 / IL-2Rα) — 7-Skill End-to-End

IL-2 cytokine (121aa) interface with IL-2Rα. Full pipeline with all 7 skills. Exposed a showstopper.

**What happened:**

- Conservation analysis: WT already at ESM-2 consensus (0 suggested mutations)
- Tight constraints: 2 redesign positions → 1 unique design (F42Y)
- Stability (interface mode): mean_LLR = -3.33 → **PASS** (relaxed gate)
- Fold-back: WT pLDDT ≈ 45, TM ≈ 0.15... wait, what?

**The showstopper:**

ESMFold couldn't fold the wildtype IL-2. pLDDT of 45 across the entire chain. The agent initially reported TM ≈ 0.15, which triggered a debug investigation.

**The debug:**

- tmtools re-alignment gave TM ≈ 0.49 (the original 0.15 was a comparison methodology error)
- But pLDDT ≈ 45 is still genuinely poor — ESMFold doesn't form the disulfide bonds (IL-2 has 3 cysteines)
- Secondary structure is correct (76.9% helix, 4-helix bundle topology matches)
- The protein folds approximately right but the confidence is low

**What was built:**

- prot-structure v4: WT sanity check — before running fold-back on designs, first check if ESMFold can fold the wildtype
- Three gate modes: HARD_GATE (WT pLDDT ≥ 60, TM ≥ 0.7), DIAGNOSTIC (report but don't reject), SKIP (ESMFold unusable)
- prot-developability v2: WT-relative mode — only flag issues that the design introduces beyond what the WT already has

**Secondary finding:**

IL-2's inherent developability flags (odd cysteines, instability index > 40) would reject every design including the wildtype under absolute thresholds. The WT-relative mode fixed this: shared flags are downgraded to INFO, only NEW flags are blocking.

## The Pipeline

```
prot-structure (quality check)
│
prot-interface (hotspots + burial)
prot-msa (conservation)
│
│  ├──── constraints ─────────────────────┘
│
prot-seqdesign (LigandMPNN + ESM-IF)
│
prot-stability (ESM-2 LLR, interface mode + conservation calibration)
│
prot-structure (fold-back, WT sanity check)
│
prot-interface (contact preservation)
│
prot-developability (red-light table, WT-relative)
│
prot-dbtl (composite scoring → decision → next cycle)
```

## Practical Lessons

### For agent trainers

**Don't teach what the agent already knows.** PE scored 9/10 on methodology. We spent zero time on protein design theory and 100% on operational tool competence.

**Every skill iteration should be driven by a real failure.** prot-structure went through 4 versions. Each version was triggered by a practice run that broke something specific. Hypothetical improvements are worth less than one real bug.

**Practice runs should get progressively harder.** 13aa peptide → 89aa standard protein → 121aa disulfide-rich cytokine. Each target was chosen to stress-test something the previous target didn't cover.

**Regression tests lock in your fixes.** After 4 versions of prot-structure, it's easy to accidentally revert an earlier fix. The regression suite catches this in seconds.

### For AI + protein engineering

**Tool disagreement is more informative than tool agreement.** The 1BRS case — fold-back PASS + interface PASS + stability FAIL — taught us more about ESM-2's limitations than any number of successful runs.

**ESM-2 is an evolutionary prior, not a stability oracle.** It captures what evolution has explored, not what's thermodynamically possible. Interface positions are conserved, so interface mutations always score poorly. Use conservation context (prot-msa) and burial analysis (prot-interface) to distinguish "ESM-2 bias" from "genuinely destabilizing."

**Always check the wildtype first.** If your evaluation tool can't handle the wildtype, it can't evaluate your designs. This applies to fold-back (ESMFold), developability (inherent protein quirks), and stability (WT LLR baseline).

**The constraint problem is real.** When >50% of interface positions are hotspots, there's nowhere left to design. All three practice runs hit this: 10/11 positions frozen (1YCR), 17/19 frozen (1BRS), 20/22 frozen (1Z92). Backbone redesign or relaxed constraint strategies are needed for tight interfaces.

## Known Limitations

- **No physics-based ΔΔG.** FoldX and Rosetta require academic licenses. ESM-2 LLR is the only stability proxy, mitigated by conservation calibration and burial context.
- **Single-chain fold-back only.** ESMFold predicts monomers. Complex structure prediction (AF2-Multimer) is not integrated.
- **ESMFold weak on disulfide-rich proteins.** Mitigated by WT sanity check, not solved.
- **Tight interfaces limit diversity.** When most positions are hotspots, LigandMPNN produces very few unique designs. Needs backbone-level redesign tools (e.g., RFdiffusion) in future work.

## Repository Structure

```
├── README.md                         # This file
├── training-summary.md               # Complete training narrative
├── skills/
│   ├── prot-seqdesign/SKILL.md       # LigandMPNN + ESM-IF (v2)
│   ├── prot-structure/SKILL.md       # ESMFold + fold-back + WT sanity (v4)
│   ├── prot-stability/SKILL.md       # ESM-2 LLR + interface mode (v2)
│   ├── prot-interface/SKILL.md       # Hotspots + burial + constraints (v2)
│   ├── prot-msa/SKILL.md             # Conservation + LLR calibration (v1)
│   ├── prot-developability/SKILL.md  # Red-light screening + WT-relative (v2)
│   └── prot-dbtl/SKILL.md            # DBTL cycle orchestration (v1)
├── tests/regression/
│   ├── run_all.py
│   ├── test_1YCR_short_peptide.py
│   ├── test_1BRS_interface_stability_conflict.py
│   └── test_1Z92_wt_sanity.py
├── docs/
│   ├── PLAYBOOK.en.md                # Operational conventions
│   ├── skill-review.en.md            # Validated capabilities summary
│   └── practice/                     # Practice run reports
│       ├── 1YCR/                     # Short peptide edge cases
│       ├── 1BRS/                     # Tool disagreement diagnosis
│       └── 1Z92/                     # 7-skill end-to-end + WT sanity debug
└── LICENSE                           # Apache 2.0
```

## Technical Stack

- **Platform:** OpenClaw
- **Model:** Claude (Anthropic)
- **Infrastructure:** NixOS on WSL2, Docker with GPU passthrough
- **Protein tools:** LigandMPNN, ESMFold (API), ESM-2 (650M), ESM-IF, fair-esm, BioPython, tmtools
- **Related project:** [ChemicalExpert training](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training) (drug discovery, same methodology)

## License

Apache 2.0. See `LICENSE`.

## Citation

```bibtex
@misc{proteinengineer2026,
  title={Teaching AI Agents to Engineer Proteins: A Systematic Training Methodology},
  author={Heng Gao},
  year={2026},
  url={https://github.com/hg125chinese-sketch/openclaw-proteinengineer-training}
}
```
