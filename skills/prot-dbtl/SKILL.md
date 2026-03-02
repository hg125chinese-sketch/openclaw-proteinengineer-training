---
name: prot-dbtl
description: Design-Build-Test-Learn cycle orchestration for protein engineering. Coordinates prot-seqdesign, prot-structure, prot-stability, and prot-interface into structured iterative cycles. Covers cycle planning, multi-skill pipeline execution, candidate selection with composite scoring, experimental feedback integration, and constraint refinement between cycles. Triggers whenever you need to run a complete design cycle, plan a multi-step protein optimization campaign, or integrate experimental results into the next round of computational design.
homepage: https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training
metadata: { "openclaw": { "emoji": "🔄", "requires": { "bins": ["python3"], "python": ["numpy", "pandas", "biopython", "torch", "fair-esm", "requests"] } } }
---

# Protein DBTL: Design-Build-Test-Learn Cycle Orchestration

Individual skills generate designs, predict structures, score stability, and analyze interfaces. This skill turns those capabilities into a repeatable, auditable optimization loop. Each cycle produces candidates with clear rationale, and experimental results feed back into the next cycle's constraints.

## When to Use

- Running a **complete design→evaluate→select pipeline** (not just one skill)
- Planning a **multi-round optimization campaign** for a protein target
- Integrating **experimental results** (expression, binding, activity) into computational constraints
- Orchestrating **multi-skill workflows** where outputs of one skill feed into another
- Making **go/no-go decisions** on candidates with composite scoring
- Documenting **why** each decision was made for reproducibility

## Core Philosophy

1. **Plan before compute.** Every cycle starts with a plan that specifies: what are we optimizing, what constraints apply, which skills will be used, and what does success look like. Get user approval before spending GPU hours.
2. **Each cycle must produce an actionable decision.** Not just data — a clear recommendation: "advance these 3 candidates" or "change strategy because X failed." If a cycle doesn't end with a decision, something went wrong.
3. **Log everything.** Which skills were used, in what order, with what parameters, and why. This is the audit trail that makes cycles reproducible and debuggable.
4. **Experimental feedback is gold.** When real data comes back, it overrides computational predictions. A design that scored poorly in silico but expressed well is MORE informative than one that scored well everywhere — it teaches you where your models are wrong.
5. **Constraints tighten across cycles.** Cycle 1 is exploratory (broad temperature, many positions). Cycle 2 uses what you learned to narrow the search. By Cycle 3 you should be making targeted mutations with high confidence.

## Phase 1: Cycle Directory Standard

### 1.1 Directory Layout

```
runs/<TARGET>_cycle<N>/
├── plan.md                    # Cycle plan (written FIRST)
├── step1_structure/           # prot-structure outputs
│   ├── <target>_clean.pdb
│   └── quality_report.md
├── step2_interface/           # prot-interface outputs
│   ├── interface_profile.csv
│   ├── hotspot_report.md
│   └── design_constraints.json
├── step3_seqdesign/           # prot-seqdesign outputs
│   ├── mpnn_output/
│   ├── designs_parsed.csv
│   └── esmif_crossval.csv
├── step4_stability/           # prot-stability outputs
│   ├── llr_scores.csv
│   ├── consensus.csv
│   └── gates_report.md
├── step5_foldback/            # prot-structure fold-back outputs
│   ├── foldback_results.csv
│   └── *.pdb
├── step6_interface_check/     # prot-interface post-design check
│   └── interface_comparison.md
├── candidates/                # Final candidates
│   ├── top_candidates.csv     # Composite scores + rationale
│   └── sequences.fasta
├── report.md                  # Cycle summary + decision
├── decisions.md               # What we decided and why
└── feedback/                  # Experimental results (Cycle N+1 input)
    ├── expression_data.csv
    ├── binding_data.csv
    └── notes.md
```

### 1.2 Create Cycle Directory

```python
#!/opt/conda/envs/prot/bin/python
"""Initialize a DBTL cycle directory."""
import os

def init_cycle(target, cycle_num, base_dir="runs"):
    """Create standardized cycle directory structure.

    Args:
        target: target name (e.g., "1BRS", "IL6R")
        cycle_num: cycle number (1, 2, 3, ...)
        base_dir: parent directory

    Returns:
        dict with directory paths
    """
    cycle_dir = os.path.join(base_dir, f"{target}_cycle{cycle_num}")

    dirs = {
        "root": cycle_dir,
        "structure": os.path.join(cycle_dir, "step1_structure"),
        "interface": os.path.join(cycle_dir, "step2_interface"),
        "seqdesign": os.path.join(cycle_dir, "step3_seqdesign"),
        "stability": os.path.join(cycle_dir, "step4_stability"),
        "foldback": os.path.join(cycle_dir, "step5_foldback"),
        "interface_check": os.path.join(cycle_dir, "step6_interface_check"),
        "candidates": os.path.join(cycle_dir, "candidates"),
        "feedback": os.path.join(cycle_dir, "feedback"),
    }

    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    print(f"Initialized {cycle_dir}/")
    return dirs
```

## Phase 2: Cycle Planning

### 2.1 Plan Template

```markdown
# Cycle Plan: <TARGET> Cycle <N>

## Objective
What are we optimizing? (e.g., "improve barstar affinity for barnase
while maintaining fold stability")

## Starting point
- Structure: <PDB ID or path>
- Target chain: <chain ID to redesign>
- Partner chain: <chain ID of binding partner>
- Previous cycle results: <link to Cycle N-1 report, or "first cycle">

## Constraints (from previous cycles)
- Frozen positions: <list, with reason>
- Must-have mutations: <if any, from experimental hits>
- Forbidden mutations: <if any, from experimental failures>
- Temperature: <LigandMPNN temperature>
- Number of designs: <N>

## Pipeline
1. [ ] prot-structure: quality check (or skip if structure unchanged)
2. [ ] prot-interface: hotspot analysis + constraint generation
3. [ ] prot-seqdesign: interface redesign with constraints
4. [ ] prot-stability: LLR scoring (mode: standard / interface)
5. [ ] prot-structure: fold-back validation
6. [ ] prot-interface: contact preservation check

## Success criteria
- At least N candidates pass all gates
- Composite score improvement over Cycle N-1 top candidates (if applicable)

## Approval
- [ ] User approved plan before execution
```

### 2.2 Write Plan Programmatically

```python
#!/opt/conda/envs/prot/bin/python
"""Generate a cycle plan from parameters."""

def write_cycle_plan(cycle_dir, target, cycle_num,
                      objective, pdb_path, design_chain, partner_chain,
                      frozen_positions=None, temperature=0.2,
                      n_designs=50, stability_mode="interface",
                      previous_cycle=None):
    """Write a plan.md for the cycle.

    Args:
        cycle_dir: path to cycle directory
        target: target name
        cycle_num: cycle number
        objective: what are we optimizing
        pdb_path: input structure
        design_chain: chain to redesign
        partner_chain: binding partner chain
        frozen_positions: list of positions to freeze (from previous cycles)
        temperature: LigandMPNN temperature
        n_designs: number of designs to generate
        stability_mode: "standard" or "interface"
        previous_cycle: path to previous cycle report (or None)

    Returns:
        path to plan.md
    """
    frozen_str = ", ".join(frozen_positions) if frozen_positions else "None (first cycle)"
    prev_str = previous_cycle if previous_cycle else "First cycle"

    plan = f"""# Cycle Plan: {target} Cycle {cycle_num}

## Objective
{objective}

## Starting point
- Structure: {pdb_path}
- Target chain: {design_chain}
- Partner chain: {partner_chain}
- Previous cycle results: {prev_str}

## Constraints
- Frozen positions: {frozen_str}
- Temperature: {temperature}
- Number of designs: {n_designs}
- Stability gate mode: {stability_mode}

## Pipeline
1. [ ] prot-structure: quality check
2. [ ] prot-interface: hotspot analysis + constraint generation
3. [ ] prot-seqdesign: interface redesign with constraints
4. [ ] prot-stability: LLR scoring (mode={stability_mode})
5. [ ] prot-structure: fold-back validation
6. [ ] prot-interface: contact preservation check

## Success criteria
- At least 3 candidates pass all gates
- If Cycle > 1: composite score improvement over previous cycle

## Approval
- [ ] User approved plan before execution
"""

    plan_path = os.path.join(cycle_dir, "plan.md")
    with open(plan_path, "w") as f:
        f.write(plan)

    print(f"Plan written → {plan_path}")
    return plan_path
```

## Phase 3: Pipeline Execution

### 3.1 Execution Order and Data Flow

```
plan.md (approved)
    │
    ▼
Step 1: prot-structure (quality check)
    │ → PDB + quality report
    ▼
Step 2: prot-interface (analysis)
    │ → hotspot profile + burial annotations + design constraints
    ▼
Step 3: prot-seqdesign (design)
    │ → N designed sequences + MPNN scores + ESM-IF cross-validation
    ▼
Step 4: prot-stability (scoring)
    │ → LLR scores + consensus + hard gates (mode from plan)
    │ → burial_info from Step 2 feeds into interface mode
    ▼
Step 5: prot-structure (fold-back)
    │ → TM-score + pLDDT for designs that passed stability
    ▼
Step 6: prot-interface (contact check)
    │ → interface preservation assessment
    ▼
Composite scoring + candidate selection
    │ → top_candidates.csv + sequences.fasta
    ▼
report.md + decisions.md
```

### 3.2 Composite Scoring

```python
#!/opt/conda/envs/prot/bin/python
"""Composite scoring for candidate selection."""
import pandas as pd
import numpy as np

def composite_score(candidates_df, task="interface_redesign"):
    """Compute composite score from multi-skill metrics.

    Weights depend on the task type.

    Args:
        candidates_df: DataFrame with columns from all skills
            Required: mpnn_score, esmif_ll, mean_llr, tm_score, plddt_foldback
            Optional: interface_preserved (bool), n_mutations, consensus_pct
        task: "interface_redesign", "stability_optimization", or "de_novo"

    Returns:
        DataFrame with composite_score column added, sorted best-first
    """
    df = candidates_df.copy()

    # Normalize each metric to [0, 1] (higher = better)
    def normalize(series, higher_is_better=True):
        s = series.fillna(series.median())
        if s.max() == s.min():
            return pd.Series(0.5, index=s.index)
        normed = (s - s.min()) / (s.max() - s.min())
        return normed if higher_is_better else (1 - normed)

    # Task-specific weights
    if task == "interface_redesign":
        weights = {
            "mpnn_score": 0.15,      # design confidence
            "esmif_ll": 0.15,        # cross-validation
            "mean_llr": 0.20,        # stability proxy
            "tm_score": 0.20,        # fold-back quality
            "plddt_foldback": 0.10,  # structural confidence
            "consensus_pct": 0.10,   # evolutionary compatibility
            "n_mutations": 0.10,     # fewer = safer (inverted)
        }
    elif task == "stability_optimization":
        weights = {
            "mpnn_score": 0.10,
            "esmif_ll": 0.10,
            "mean_llr": 0.35,        # stability is the objective
            "tm_score": 0.20,
            "plddt_foldback": 0.15,
            "consensus_pct": 0.05,
            "n_mutations": 0.05,
        }
    else:  # de_novo
        weights = {
            "mpnn_score": 0.20,
            "esmif_ll": 0.20,
            "mean_llr": 0.15,
            "tm_score": 0.25,
            "plddt_foldback": 0.15,
            "consensus_pct": 0.05,
            "n_mutations": 0.00,
        }

    # Compute normalized scores
    score = pd.Series(0.0, index=df.index)
    for metric, weight in weights.items():
        if metric in df.columns and weight > 0:
            higher_better = metric != "n_mutations"
            score += weight * normalize(df[metric], higher_is_better=higher_better)

    df["composite_score"] = score
    df = df.sort_values("composite_score", ascending=False)

    return df
```

### 3.3 Gate Summary

```python
#!/opt/conda/envs/prot/bin/python
"""Summarize gate results across all skills."""

def gate_summary(candidates_df):
    """Check which candidates pass all gates.

    Expected gate columns (True/False):
    - stability_pass: from prot-stability hard gates
    - foldback_pass: from prot-structure fold-back (TM >= threshold)
    - interface_pass: from prot-interface hotspot gates

    Returns:
        all_pass_df, partial_pass_df, summary_dict
    """
    gate_cols = [c for c in candidates_df.columns if c.endswith("_pass")]

    df = candidates_df.copy()
    df["all_gates_pass"] = df[gate_cols].all(axis=1)
    df["n_gates_pass"] = df[gate_cols].sum(axis=1)

    all_pass = df[df["all_gates_pass"]]
    partial = df[~df["all_gates_pass"]]

    summary = {
        "n_candidates": len(df),
        "n_all_pass": len(all_pass),
        "n_partial": len(partial),
        "gate_pass_rates": {col: df[col].mean() for col in gate_cols},
    }

    print(f"Gate summary: {summary['n_all_pass']}/{summary['n_candidates']} pass all gates")
    for col, rate in summary["gate_pass_rates"].items():
        print(f"  {col}: {rate:.0%}")

    return all_pass, partial, summary
```

## Phase 4: Experimental Feedback Integration

### 4.1 Load and Parse Experimental Results

```python
#!/opt/conda/envs/prot/bin/python
"""Parse experimental results and extract learnings for next cycle."""
import pandas as pd

def load_experimental_feedback(feedback_dir):
    """Load experimental data from the feedback directory.

    Expected files (all optional):
    - expression_data.csv: columns [name, expressed (bool), yield_mg_L (float)]
    - binding_data.csv: columns [name, kd_nM (float), fold_change (float)]
    - notes.md: free-text observations

    Returns:
        dict with parsed data
    """
    import os

    result = {}

    expr_path = os.path.join(feedback_dir, "expression_data.csv")
    if os.path.exists(expr_path):
        result["expression"] = pd.read_csv(expr_path)
        n_expr = result["expression"]["expressed"].sum()
        n_total = len(result["expression"])
        print(f"Expression: {n_expr}/{n_total} expressed")

    bind_path = os.path.join(feedback_dir, "binding_data.csv")
    if os.path.exists(bind_path):
        result["binding"] = pd.read_csv(bind_path)
        print(f"Binding data: {len(result['binding'])} entries")

    notes_path = os.path.join(feedback_dir, "notes.md")
    if os.path.exists(notes_path):
        with open(notes_path) as f:
            result["notes"] = f.read()

    return result


def extract_learnings(feedback, computational_predictions):
    """Compare experimental results with computational predictions.

    Identifies:
    - True positives: predicted good AND experimentally good
    - False positives: predicted good BUT experimentally bad
    - False negatives: predicted bad BUT experimentally good (most informative!)
    - True negatives: predicted bad AND experimentally bad

    Args:
        feedback: from load_experimental_feedback()
        computational_predictions: DataFrame with columns [name, composite_score, ...]

    Returns:
        learnings dict with categorized results and constraint suggestions
    """
    learnings = {
        "true_positives": [],
        "false_positives": [],
        "false_negatives": [],
        "true_negatives": [],
        "constraint_suggestions": [],
    }

    # Merge expression data if available
    if "expression" in feedback:
        expr = feedback["expression"]
        merged = computational_predictions.merge(expr, on="name", how="inner")

        for _, row in merged.iterrows():
            predicted_good = row.get("composite_score", 0) > 0.5
            expressed = row.get("expressed", False)

            if predicted_good and expressed:
                learnings["true_positives"].append(row["name"])
            elif predicted_good and not expressed:
                learnings["false_positives"].append(row["name"])
                # Analyze why: which mutations are problematic?
                learnings["constraint_suggestions"].append(
                    f"FREEZE mutations in {row['name']} — predicted good but failed expression"
                )
            elif not predicted_good and expressed:
                learnings["false_negatives"].append(row["name"])
                learnings["constraint_suggestions"].append(
                    f"REVIEW gates that rejected {row['name']} — it actually expressed"
                )
            else:
                learnings["true_negatives"].append(row["name"])

    # Summary
    tp = len(learnings["true_positives"])
    fp = len(learnings["false_positives"])
    fn = len(learnings["false_negatives"])
    tn = len(learnings["true_negatives"])
    total = tp + fp + fn + tn

    if total > 0:
        learnings["accuracy"] = (tp + tn) / total
        learnings["precision"] = tp / (tp + fp) if (tp + fp) > 0 else 0
        learnings["recall"] = tp / (tp + fn) if (tp + fn) > 0 else 0
        print(f"Prediction accuracy: {learnings['accuracy']:.0%} "
              f"(precision={learnings['precision']:.0%}, recall={learnings['recall']:.0%})")
        if fn > 0:
            print(f"⚠️  {fn} false negatives — review gate thresholds!")
        if fp > 0:
            print(f"⚠️  {fp} false positives — review scoring model!")

    return learnings
```

### 4.2 Generate Next-Cycle Constraints

```python
#!/opt/conda/envs/prot/bin/python
"""Generate updated constraints for the next cycle based on learnings."""

def next_cycle_constraints(learnings, current_constraints, hotspot_profile=None):
    """Update design constraints based on experimental feedback.

    Rules:
    - False positives → tighten gates or freeze their mutation positions
    - False negatives → relax gates that rejected them
    - Consistent failures at a position → freeze that position
    - Consistent successes at a position → keep it in redesign set

    Args:
        learnings: from extract_learnings()
        current_constraints: dict with current frozen/redesign lists
        hotspot_profile: from prot-interface (optional, for position mapping)

    Returns:
        updated constraints dict for next cycle
    """
    new_constraints = {
        "frozen": list(current_constraints.get("frozen", [])),
        "redesign": list(current_constraints.get("redesign", [])),
        "temperature_adjustment": None,
        "gate_mode_suggestion": current_constraints.get("gate_mode", "interface"),
        "notes": [],
    }

    # From false positives: freeze positions that failed experimentally
    if learnings.get("false_positives"):
        new_constraints["notes"].append(
            f"False positives detected: {learnings['false_positives']}. "
            f"Consider tightening stability gates or freezing their mutation positions."
        )
        # Suggest lowering temperature for more conservative designs
        new_constraints["temperature_adjustment"] = "decrease"

    # From false negatives: relax gates
    if learnings.get("false_negatives"):
        new_constraints["notes"].append(
            f"False negatives detected: {learnings['false_negatives']}. "
            f"Review which gates rejected them and consider relaxing thresholds."
        )

    # High false positive rate → tighten everything
    if learnings.get("precision", 1) < 0.5:
        new_constraints["notes"].append(
            "Low precision (<50%). Recommend: lower temperature, "
            "switch to standard stability mode, or add physics-based scoring."
        )

    # High false negative rate → loosen gates
    if learnings.get("recall", 1) < 0.5:
        new_constraints["notes"].append(
            "Low recall (<50%). Recommend: relax LLR gates, "
            "use interface mode if not already, or increase temperature."
        )

    return new_constraints
```

## Phase 5: Decision Documentation

### 5.1 Write Cycle Report

```python
#!/opt/conda/envs/prot/bin/python
"""Generate cycle report and decision document."""
import os

def write_cycle_report(cycle_dir, target, cycle_num,
                        gate_summary_dict, top_candidates_df,
                        learnings=None, decision=None):
    """Write the cycle summary report and decisions.

    Args:
        cycle_dir: cycle directory path
        target: target name
        cycle_num: cycle number
        gate_summary_dict: from gate_summary()
        top_candidates_df: final ranked candidates
        learnings: from extract_learnings() (if experimental data available)
        decision: string describing the decision

    Returns:
        paths to report.md and decisions.md
    """
    # Report
    report = f"# {target} Cycle {cycle_num} — Report\n\n"
    report += f"## Gate Summary\n"
    report += f"- Candidates evaluated: {gate_summary_dict['n_candidates']}\n"
    report += f"- All gates pass: {gate_summary_dict['n_all_pass']}\n"
    for gate, rate in gate_summary_dict.get("gate_pass_rates", {}).items():
        report += f"- {gate}: {rate:.0%}\n"

    report += f"\n## Top Candidates\n"
    if len(top_candidates_df) > 0:
        cols = [c for c in ["name", "composite_score", "mean_llr", "tm_score",
                           "plddt_foldback", "n_mutations"] if c in top_candidates_df.columns]
        report += top_candidates_df[cols].head(10).to_markdown(index=False)
    else:
        report += "No candidates passed all gates.\n"

    if learnings:
        report += f"\n\n## Experimental Feedback Analysis\n"
        report += f"- Accuracy: {learnings.get('accuracy', 'N/A')}\n"
        report += f"- True positives: {learnings.get('true_positives', [])}\n"
        report += f"- False positives: {learnings.get('false_positives', [])}\n"
        report += f"- False negatives: {learnings.get('false_negatives', [])}\n"

    report_path = os.path.join(cycle_dir, "report.md")
    with open(report_path, "w") as f:
        f.write(report)

    # Decisions
    if decision is None:
        if gate_summary_dict["n_all_pass"] >= 3:
            decision = (f"Advance top 3 candidates to experimental validation. "
                       f"Cycle {cycle_num} successful.")
        elif gate_summary_dict["n_all_pass"] > 0:
            decision = (f"Only {gate_summary_dict['n_all_pass']} candidates pass. "
                       f"Advance them but consider adjusting strategy for Cycle {cycle_num+1}.")
        else:
            decision = (f"No candidates pass all gates. "
                       f"Review pipeline: check if constraints are too tight, "
                       f"temperature too low, or structure quality insufficient. "
                       f"Recommended: adjust constraints and re-run as Cycle {cycle_num+1}.")

    decisions_md = f"# {target} Cycle {cycle_num} — Decisions\n\n"
    decisions_md += f"## Decision\n{decision}\n\n"
    decisions_md += f"## Skills used\n"
    decisions_md += "1. prot-structure (quality check + fold-back)\n"
    decisions_md += "2. prot-interface (hotspot + constraints + contact check)\n"
    decisions_md += "3. prot-seqdesign (LigandMPNN + ESM-IF)\n"
    decisions_md += "4. prot-stability (ESM-2 LLR + consensus)\n"

    decisions_path = os.path.join(cycle_dir, "decisions.md")
    with open(decisions_path, "w") as f:
        f.write(decisions_md)

    print(f"Report → {report_path}")
    print(f"Decisions → {decisions_path}")

    return report_path, decisions_path
```

## Phase 6: Multi-Cycle Strategy

### 6.1 Cycle Progression Pattern

```
Cycle 1: EXPLORATION
├── Temperature: 0.2-0.3 (moderate diversity)
├── Redesign: all tolerant + caution positions
├── Stability mode: interface
├── Goal: find promising mutation patterns
└── Expected: 3-10 candidates pass gates

Cycle 2: REFINEMENT (uses Cycle 1 feedback)
├── Temperature: 0.1-0.2 (lower, more focused)
├── Redesign: only positions where Cycle 1 hits occurred
├── Freeze: positions where Cycle 1 mutations failed
├── Stability mode: interface (with burial)
├── Goal: optimize promising patterns
└── Expected: 5-15 candidates pass gates (narrower but better)

Cycle 3: PRECISION (uses Cycle 2 feedback)
├── Temperature: 0.1 (very focused)
├── Redesign: 2-3 specific positions
├── Combinations: try different mutations at validated positions
├── Stability mode: standard (for final candidates)
├── Goal: identify best 1-3 candidates for scale-up
└── Expected: high-confidence final candidates
```

### 6.2 When to Change Strategy

```
After each cycle, ask:
│
├── Did ≥3 candidates pass all gates?
│   ├── YES → Continue to next cycle with tightened constraints
│   └── NO →
│       ├── Which gate blocked the most?
│       │   ├── Stability → switch to interface mode / relax LLR
│       │   ├── Fold-back → structure quality issue → try different backbone
│       │   └── Interface → too many hotspots → expand redesign set
│       │
│       ├── Is the interface too constrained (>80% hotspots)?
│       │   ├── YES → Consider designing away from interface
│       │   │         or use backbone redesign (prot-backbone-gen, future)
│       │   └── NO → Adjust temperature or constraints
│       │
│       └── Did previous cycle's experiments contradict predictions?
│           ├── YES → Recalibrate: which metric was wrong?
│           └── NO → Tighten constraints based on consistent patterns
```

## Phase 7: VRAM Budget & Timing (Full Cycle)

| Step | VRAM peak | Wall time | Notes |
|------|-----------|-----------|-------|
| Structure quality check | 0 GB | ~5 sec | CPU only |
| Interface analysis + hotspot | ~2 GB | ~15 sec | ESM-2 alanine scan |
| LigandMPNN (50 designs) | ~2 GB | ~30 sec | GPU |
| ESM-IF cross-validation (Top 10) | ~2 GB | ~2 min | GPU |
| Stability scoring (Top 10) | ~2 GB | ~10 sec | ESM-2 LLR |
| Fold-back (5 designs via API) | 0 GB | ~1-2 min | ESMFold API |
| Interface contact check | 0 GB | ~5 sec | CPU |
| **Total per cycle** | **~2 GB peak** | **~5-10 min** | Sequential steps |

## Failure Modes

1. **Skipping the plan.** The plan forces you to think about constraints and success criteria before burning compute. Without it, you'll redesign the wrong positions or use the wrong gates.

2. **Not logging which skills were used.** When a cycle fails, you need to know what happened. decisions.md exists for this. Don't skip it.

3. **Ignoring false negatives.** When experimental data shows a design worked that your pipeline rejected, that's the most valuable signal — it tells you where your models are wrong. Don't dismiss it.

4. **Over-tightening after failures.** If Cycle 1 produces 0 candidates, the instinct is to make everything stricter. Usually the opposite is needed: relax gates, increase temperature, or expand the redesign set.

5. **Running too many cycles without experimental data.** Computational cycles are cheap but they can drift. Every 2-3 computational cycles should include at least one round of experimental validation.

6. **Changing multiple variables between cycles.** If you change temperature AND constraints AND stability mode between cycles, you can't attribute improvements. Change one thing at a time when possible.

## One-Sentence Rule

**Plan every cycle before executing, log every decision with rationale, and when experimental data comes back — let it override your computational predictions and reshape your constraints.**
