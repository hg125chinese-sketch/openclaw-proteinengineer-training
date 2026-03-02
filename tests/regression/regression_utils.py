import json
import os
import re
from dataclasses import dataclass


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


def load_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def assert_regex(text: str, pattern: str, msg: str):
    if re.search(pattern, text, flags=re.MULTILINE) is None:
        raise AssertionError(f"{msg} | missing pattern: {pattern}")


def assert_contains(text: str, needle: str, msg: str):
    if needle not in text:
        raise AssertionError(f"{msg} | missing substring: {needle}")


def run_test(test_name: str, fn, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    result = {
        "test": test_name,
        "passed": False,
        "checks": [],
        "error": None,
    }
    try:
        checks = fn()
        # checks: list[Check]
        passed = all(c.ok for c in checks)
        result["passed"] = passed
        result["checks"] = [c.__dict__ for c in checks]
    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}"

    # write individual result
    out_path = os.path.join(out_dir, f"{test_name}.result.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    # print human-friendly
    if result["passed"]:
        print(f"PASS {test_name}")
    else:
        print(f"FAIL {test_name}")
        if result["error"]:
            print("  error:", result["error"])
        for c in result.get("checks", []):
            if not c.get("ok"):
                print(f"  - {c.get('name')}: {c.get('detail')}")

    return result
