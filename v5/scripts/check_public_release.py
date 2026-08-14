"""Check the curated V5 repository boundary without research inputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {
    ".rda", ".rdata", ".rds", ".h5", ".h5ad", ".loom", ".bam",
    ".fastq", ".fq", ".xlsx", ".zip",
}
FORBIDDEN_HEADERS = {
    "sample_accession", "sample_id", "subject_id", "cell_id", "donor_id",
    "geo_accession", "count_column",
}
FORBIDDEN_MARKERS = (
    "github_pat_", "ghp_", "begin rsa private key",
    "begin openssh private key", "begin private key",
)


def files() -> list[Path]:
    return [
        path for path in ROOT.rglob("*")
        if path.is_file() and ".git" not in path.parts and path.suffix.lower() != ".pyc"
    ]


def main() -> int:
    errors: list[str] = []
    release_files = files()
    for path in release_files:
        relative = path.relative_to(ROOT)
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"forbidden file type: {relative}")
        if path.stat().st_size >= 100 * 1024 * 1024:
            errors.append(f"file exceeds GitHub single-file limit: {relative}")
        if path.name != "check_public_release.py" and path.suffix.lower() in {".md", ".cff", ".json", ".py", ".txt"}:
            text = path.read_text(encoding="utf-8", errors="replace").lower()
            if any(marker in text for marker in FORBIDDEN_MARKERS):
                errors.append(f"credential marker in {relative}")

    for path in release_files:
        if path.suffix.lower() != ".csv":
            continue
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            header = {column.strip().lower() for column in next(csv.reader(handle), [])}
        leaked = sorted(header & FORBIDDEN_HEADERS)
        if leaked:
            errors.append(f"row-level header in {path.relative_to(ROOT)}: {', '.join(leaked)}")

    report_path = ROOT / "99_QA" / "qa_report.json"
    if not report_path.is_file():
        errors.append("missing 99_QA/qa_report.json")
    else:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        if report.get("status") != "PASS" or report.get("failure_count") != 0:
            errors.append("V5 QA report is not PASS with zero failures")

    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"OK: curated V5 public release checks passed for {len(release_files)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
