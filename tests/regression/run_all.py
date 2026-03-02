#!/opt/conda/envs/prot/bin/python
"""Run all regression tests and write summary.json.

Each test writes <test>.result.json. This script aggregates into summary.json.
"""

import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))

TESTS = [
    "test_1YCR_short_peptide.py",
    "test_1BRS_interface_stability_conflict.py",
    "test_1Z92_wt_sanity.py",
]


def main():
    results = []
    for t in TESTS:
        path = os.path.join(HERE, t)
        p = subprocess.run(["/opt/conda/envs/prot/bin/python", path], capture_output=True, text=True)
        print(p.stdout, end="")
        if p.returncode != 0:
            print(p.stderr)

        # load per-test result
        name = t.replace(".py", "")
        res_path = os.path.join(HERE, f"{name}.result.json")
        if os.path.exists(res_path):
            results.append(json.load(open(res_path)))
        else:
            results.append({"test": name, "passed": False, "error": "no result json produced"})

    summary = {
        "n_tests": len(results),
        "n_pass": sum(1 for r in results if r.get("passed")),
        "n_fail": sum(1 for r in results if not r.get("passed")),
        "results": results,
    }

    out_path = os.path.join(HERE, "summary.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nWrote {out_path}")
    if summary["n_fail"] > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
