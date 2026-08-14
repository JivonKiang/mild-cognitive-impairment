"""Audit GSE63063 as a complete GSE63060/GSE63061 GSM aggregate."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_public_discovery_20260812\gds_mci_blood_summary.json"
)
DEFAULT_RESULTS_OUTPUT = ROOT / "03_RESULTS" / "gse63063_overlap_audit.json"
DEFAULT_SOURCE_OUTPUT = ROOT / "04_FIGURES" / "source_data" / "gse63063_overlap_audit.json"
DEVELOPMENT_TARGETS = ("GSE63060", "GSE63061")
TARGETS = (*DEVELOPMENT_TARGETS, "GSE63063", "GSE140829")


def sha256(path: Path) -> str:
    """Return the SHA-256 digest for a file without loading it all into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_source_provenance(summary_path: Path) -> dict[str, str]:
    """Return filesystem provenance for the validated GEO summary cache."""
    modified_utc = datetime.fromtimestamp(
        summary_path.stat().st_mtime,
        tz=timezone.utc,
    ).isoformat().replace("+00:00", "Z")
    return {
        "metadata_source": str(summary_path),
        "metadata_source_sha256": sha256(summary_path),
        "metadata_source_modified_utc": modified_utc,
    }


def load_and_validate_series(summary_path: Path) -> dict[str, dict[str, Any]]:
    """Load target accession metadata and validate its cache schema and GSM sets."""
    if not summary_path.is_file():
        raise FileNotFoundError(f"official GEO summary cache not found: {summary_path}")
    try:
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"official GEO summary cache is not valid JSON: {summary_path}") from error

    result = payload.get("result") if isinstance(payload, dict) else None
    if not isinstance(result, dict):
        raise ValueError("official GEO summary cache must contain an object at 'result'")

    series: dict[str, dict[str, Any]] = {}
    for record_key, record in result.items():
        if not isinstance(record, dict):
            continue
        accession = record.get("accession")
        if accession in TARGETS:
            if not isinstance(accession, str):
                raise ValueError(f"cache record {record_key!r} has an invalid accession")
            if accession in series:
                raise ValueError(f"duplicate accession record in cache: {accession}")
            series[accession] = record

    missing_targets = sorted(set(TARGETS) - set(series))
    if missing_targets:
        raise ValueError(f"official GEO summary cache is missing target accessions: {', '.join(missing_targets)}")

    for accession, record in series.items():
        samples = record.get("samples")
        if not isinstance(samples, list):
            raise ValueError(f"{accession} cache record must contain a list at 'samples'")
        if not samples:
            raise ValueError(f"{accession} cache record has an empty sample set")
        sample_accessions: set[str] = set()
        for index, sample in enumerate(samples):
            if not isinstance(sample, dict):
                raise ValueError(f"{accession} sample {index} must be an object")
            sample_accession = sample.get("accession")
            if not isinstance(sample_accession, str) or not sample_accession.startswith("GSM") or not sample_accession[3:].isdigit():
                raise ValueError(f"{accession} sample {index} has a missing or invalid GSM accession: {sample_accession!r}")
            if sample_accession in sample_accessions:
                raise ValueError(f"{accession} cache record has duplicate GSM accession: {sample_accession}")
            sample_accessions.add(sample_accession)
    return series


def build_audit(series: dict[str, dict[str, Any]], provenance: Mapping[str, str]) -> dict[str, Any]:
    """Build the accession-only independence audit from validated GEO metadata."""
    missing_targets = sorted(set(TARGETS) - set(series))
    if missing_targets:
        raise ValueError(f"validated series is missing target accessions: {', '.join(missing_targets)}")
    gsm = {
        accession: {sample["accession"] for sample in series[accession]["samples"]}
        for accession in TARGETS
    }
    development_union = gsm["GSE63060"] | gsm["GSE63061"]
    aggregate = gsm["GSE63063"]
    gse140829_shared = aggregate & gsm["GSE140829"]
    overlap = []
    for accession in DEVELOPMENT_TARGETS:
        shared = gsm[accession] & aggregate
        overlap.append(
            {
                "development_series": accession,
                "development_sample_count": len(gsm[accession]),
                "shared_gsm_count": len(shared),
                "fraction_of_development_series": len(shared) / len(gsm[accession]),
                "fraction_of_gse63063": len(shared) / len(aggregate),
                "shared_gsm_examples": sorted(shared)[:10],
            }
        )

    expected_shared_counts = {"GSE63060": 329, "GSE63061": 388}
    observed_shared_counts = {
        row["development_series"]: row["shared_gsm_count"] for row in overlap
    }
    observed_fractions = {
        row["development_series"]: row["fraction_of_development_series"] for row in overlap
    }
    if (
        observed_shared_counts != expected_shared_counts
        or observed_fractions != {accession: 1.0 for accession in DEVELOPMENT_TARGETS}
        or len(development_union) != 717
        or len(aggregate - development_union) != 0
        or len(development_union - aggregate) != 0
        or len(gse140829_shared) != 0
    ):
        raise ValueError("GSE63063 aggregate invariants do not match the audited boundary")

    return {
        "status": "AUDITED_COMPLETE_DEVELOPMENT_AGGREGATE",
        "decision": "EXCLUDE_FROM_INDEPENDENT_VALIDATION_AND_MODELING",
        "decision_rationale": (
            "The official GEO summary cache shows that every GSE63060 and GSE63061 GSM record "
            "is contained in GSE63063, with no additional GSE63063 GSM records. GSE63063 is a "
            "complete aggregate of the AddNeuroMed development series, not a new cohort."
        ),
        "gse63063_sample_count": len(aggregate),
        "development_union_sample_count": len(development_union),
        "gse63063_extra_gsm_count": len(aggregate - development_union),
        "development_union_missing_from_gse63063_count": len(development_union - aggregate),
        "gse140829_shared_gsm_count": len(gse140829_shared),
        "gse140829_shared_gsm_examples": sorted(gse140829_shared)[:10],
        "development_overlap": overlap,
        "series": {
            accession: {
                "title": series[accession].get("title"),
                "platform": series[accession].get("gpl"),
                "sample_count": len(gsm[accession]),
                "sample_accession_hash": hashlib.sha256(
                    "\n".join(sorted(gsm[accession])).encode("utf-8")
                ).hexdigest(),
                "pubmed_ids": series[accession].get("pubmedids", []),
            }
            for accession in TARGETS
        },
        "analysis_boundary": (
            "No GSE63063 expression values or phenotype labels are used for model fitting, feature "
            "selection, preprocessing, calibration, thresholding, adaptation, or independent validation. "
            "The audit is an accession-level sample-identity check only."
        ),
        "official_sources": {
            accession: f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
            for accession in TARGETS
        },
        "metadata_source": provenance["metadata_source"],
        "metadata_source_sha256": provenance["metadata_source_sha256"],
        "metadata_source_modified_utc": provenance["metadata_source_modified_utc"],
    }


def write_audit(audit: dict[str, Any], results_output: Path, source_output: Path) -> None:
    """Write byte-identical result and source-data copies of an audit."""
    rendered = json.dumps(audit, indent=2) + "\n"
    for destination in (results_output, source_output):
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(rendered, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY, help="official GEO summary JSON cache")
    parser.add_argument("--results-output", type=Path, default=DEFAULT_RESULTS_OUTPUT, help="audit result JSON path")
    parser.add_argument("--source-output", type=Path, default=DEFAULT_SOURCE_OUTPUT, help="source-data JSON copy path")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    audit = build_audit(load_and_validate_series(args.summary), build_source_provenance(args.summary))
    write_audit(audit, args.results_output, args.source_output)
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
