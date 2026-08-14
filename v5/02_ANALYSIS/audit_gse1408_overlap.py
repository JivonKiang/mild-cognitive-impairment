"""Audit GSE140829/830/831 overlap before treating a series as independent."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_public_discovery_20260812\gds_mci_blood_summary.json"
)
OUT = ROOT / "03_RESULTS"
SOURCE_OUT = ROOT / "04_FIGURES" / "source_data"
TARGETS = ("GSE140829", "GSE140830", "GSE140831")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_series() -> dict[str, dict[str, object]]:
    payload = json.loads(SUMMARY.read_text(encoding="utf-8"))
    result = payload.get("result", {})
    series = {
        value.get("accession"): value
        for value in result.values()
        if isinstance(value, dict) and value.get("accession") in TARGETS
    }
    if set(series) != set(TARGETS):
        raise RuntimeError(f"missing target series in official NCBI summary: {TARGETS}")
    return series


def label_counts(samples: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for sample in samples:
        title = sample["title"]
        label = title.split(",", 2)[1].strip() if "," in title else "UNPARSED"
        counts[label] = counts.get(label, 0) + 1
    return dict(sorted(counts.items()))


def main() -> int:
    if not SUMMARY.exists():
        raise FileNotFoundError(SUMMARY)
    series = load_series()
    accessions = {
        accession: {sample["accession"] for sample in value["samples"]}
        for accession, value in series.items()
    }
    overlap_rows = []
    for index, left in enumerate(TARGETS):
        for right in TARGETS[index + 1 :]:
            intersection = sorted(accessions[left] & accessions[right])
            overlap_rows.append(
                {
                    "left": left,
                    "right": right,
                    "left_sample_count": len(accessions[left]),
                    "right_sample_count": len(accessions[right]),
                    "shared_gsm_count": len(intersection),
                    "shared_fraction_of_left": len(intersection) / len(accessions[left]),
                    "shared_fraction_of_right": len(intersection) / len(accessions[right]),
                    "shared_gsm_examples": intersection[:10],
                }
            )

    records = {
        accession: {
            "title": series[accession]["title"],
            "platform": series[accession]["gpl"],
            "assay": series[accession]["gdstype"],
            "sample_count": len(accessions[accession]),
            "label_counts": label_counts(series[accession]["samples"]),
            "sample_accession_hash": hashlib.sha256(
                "\n".join(sorted(accessions[accession])).encode("utf-8")
            ).hexdigest(),
        }
        for accession in TARGETS
    }
    audit = {
        "status": "AUDITED_NONINDEPENDENT_CANDIDATES",
        "decision": {
            "GSE140829": "retained as the existing external transport assessment",
            "GSE140830": "excluded from primary MCI-versus-HC validation; no MCI samples and same parent study family",
            "GSE140831": "excluded as non-independent; contains all GSE140829 and GSE140830 GSM records",
        },
        "series": records,
        "pairwise_overlap": overlap_rows,
        "interpretation_boundary": (
            "GSE140830 and GSE140831 must not be counted as independent replication. "
            "GSE140831 is an aggregate/overlapping series, while GSE140830 has no MCI label. "
            "The overlap audit uses no expression values or phenotype labels for model fitting, thresholding, "
            "calibration, feature selection, or adaptation."
        ),
        "official_sources": {
            accession: f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
            for accession in TARGETS
        },
        "metadata_source": str(SUMMARY),
        "metadata_source_sha256": sha256(SUMMARY),
        "access_date": "2026-08-13",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    for destination in [OUT / "gse1408_overlap_audit.json", SOURCE_OUT / "gse1408_overlap_audit.json"]:
        destination.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
