# ProteinEngineer Playbook (Public, English Extract)

This is a **sanitized** extract of the internal playbook, focused on reproducible, verifiable protein-engineering workflows.

## Verifiability Standard

When asked to “check docs / provide verifiable instructions / follow the playbook”, follow this pattern:

1. **Use QMD first** (select the most relevant collection)
2. **Show evidence** in your final response:
   - the `qmd://...` reference **with line numbers**
   - a short quoted excerpt (with line numbers)
   - exact copy/paste commands
3. If no collection covers the topic, say so and propose adding a new entry.

## Environment Convention

- Use the project’s **protein** Python environment for protein work.
- Do not silently switch environments.

## Engineering Hygiene Rule

- After changing any skill document or pipeline logic, run:

```bash
./scripts/check.sh
```

Only commit if **all regression tests PASS**.
