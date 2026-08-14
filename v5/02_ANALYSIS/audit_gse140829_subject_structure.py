from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse140829_task8\GSE140829_series_matrix.fresh.txt.gz"
)
OUT = ROOT / "03_RESULTS"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean(value: str) -> str:
    return value.strip().strip('"')


def parse() -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    with gzip.open(SOURCE, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw in handle:
            if not raw.startswith("!Sample_"):
                continue
            fields = next(csv.reader([raw.rstrip("\r\n")], delimiter="\t"), [])
            if fields:
                rows.setdefault(fields[0], []).append([clean(item) for item in fields[1:]])
    return {key: values[0] for key, values in rows.items() if values}


def label(title: str) -> str:
    lower = title.casefold()
    # Only accept a diagnosis token in the comma-delimited title field.
    # Series suffixes such as ``[ad_mci]`` are not sample labels.
    match = re.search(r"(?:^|,\s*)(control|mci|ad)(?:\s*,|$)", lower)
    if match:
        return {"control": "HC", "mci": "MCI", "ad": "AD"}[match.group(1)]
    return "OTHER"


def candidate_subject_id(title: str) -> str:
    # This is intentionally a candidate only. It is never promoted without
    # an accession-to-subject map from the originating study.
    fields = [field.strip() for field in title.strip().split(",")]
    sample_field = fields[2] if len(fields) >= 3 else ""
    token = re.search(
        r"([A-Za-z0-9][A-Za-z0-9.-]{5,})(?:_[A-Za-z])?(?:\s*\[|$)",
        sample_field,
    )
    return token.group(1) if token else "UNRESOLVED"


def main() -> int:
    rows = parse()
    titles = rows.get("!Sample_title", [])
    accessions = rows.get("!Sample_geo_accession", [])
    n = min(len(titles), len(accessions))
    records = []
    for accession, title in zip(accessions[:n], titles[:n]):
        records.append(
            {
                "sample_accession": accession,
                "title": title,
                "label_from_title": label(title),
                "candidate_title_prefix": candidate_subject_id(title),
            }
        )
    labels = Counter(row["label_from_title"] for row in records)
    mci_hc = [row for row in records if row["label_from_title"] in {"MCI", "HC"}]
    prefix_counts = Counter(row["candidate_title_prefix"] for row in mci_hc)
    repeated = sum(count - 1 for count in prefix_counts.values() if count > 1)
    cross_label = {}
    for prefix in prefix_counts:
        cross_label[prefix] = sorted({row["label_from_title"] for row in records if row["candidate_title_prefix"] == prefix})
    audit = {
        "status": "AUDIT_ONLY_UNRESOLVED_SUBJECT_INDEPENDENCE",
        "source": str(SOURCE),
        "source_sha256": sha256(SOURCE),
        "series_sample_rows": len(records),
        "title_label_counts": dict(sorted(labels.items())),
        "primary_mci_hc_rows": len(mci_hc),
        "candidate_prefix_unique_mci_hc": len(prefix_counts),
        "candidate_prefix_repeated_rows_mci_hc": repeated,
        "candidate_prefixes_crossing_labels": {
            key: value for key, value in cross_label.items() if len(value) > 1
        },
        "interpretation": (
            "Candidate title prefixes are descriptive diagnostics only. "
            "The external confidence interval must not be reported as subject-level "
            "inference until the originating study provides a validated sample-to-subject map."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gse140829_subject_structure_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    with (OUT / "gse140829_subject_structure_records.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys() if records else ["sample_accession"])
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
