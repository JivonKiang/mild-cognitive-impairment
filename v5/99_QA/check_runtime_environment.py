"""Verify the direct Python dependencies used by the V5 analysis package."""

from __future__ import annotations

import json
import sys
from importlib.metadata import version


REQUIRED = {
    "matplotlib": "3.10.8",
    "numpy": "2.2.6",
    "pandas": "2.3.3",
    "Pillow": "12.2.0",
    "scikit-learn": "1.7.2",
    "scipy": "1.15.3",
    "statsmodels": "0.14.6",
}


def main() -> int:
    observed = {package: version(package) for package in REQUIRED}
    mismatches = {
        package: {"expected": expected, "observed": observed[package]}
        for package, expected in REQUIRED.items()
        if observed[package] != expected
    }
    print(
        json.dumps(
            {
                "status": "PASS" if not mismatches else "FAIL",
                "python_version": sys.version.split()[0],
                "observed": observed,
                "mismatches": mismatches,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
