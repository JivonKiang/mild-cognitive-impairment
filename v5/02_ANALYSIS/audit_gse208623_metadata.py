"""Audit GSE208623 metadata without reading methylation beta values.

The cohort is retained only as independently sourced cross-omics context. This
module intentionally cannot fit transcriptome models, analyse CpGs, or infer
mechanisms. It writes only aggregate public metadata; the separate source-locked
association workflow parses the official header internally when it needs row-level
alignment.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "01_DATA_AUDIT" / "gse208623"
HEADER_PATH = INPUT_DIR / "GSE208623_series_matrix_header.gz"
BETA_PATH = INPUT_DIR / "GSE208623_Norm_beta.csv.gz"
OUTPUT_DIR = ROOT / "03_RESULTS"
BETA_SOURCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208623/suppl/"
    "GSE208623_Norm_beta.csv.gz"
)
OFFICIAL_BETA_EXPECTED_BYTES = 363_136_110


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_metadata(path: Path) -> list[dict[str, str]]:
    fields: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *values = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            if key == "!Sample_characteristics_ch1":
                characteristics.append(values)
            else:
                fields[key] = values
    titles = fields.get("!Sample_title", [])
    gsm = fields.get("!Sample_geo_accession", [])
    if not titles or len(titles) != 60 or len(gsm) != len(titles):
        raise ValueError("GSE208623 official header does not contain 60 aligned sample records")
    if any(len(values) != len(titles) for values in characteristics):
        raise ValueError("GSE208623 characteristic row is not aligned to sample metadata")
    rows: list[dict[str, str]] = []
    for index, title in enumerate(titles):
        sex = "MISSING"
        for values in characteristics:
            name, separator, value = values[index].partition(":")
            if separator and name.strip().casefold() == "gender":
                sex = value.strip().casefold()
        title_folded = title.casefold()
        if title_folded.startswith("mci subject"):
            endpoint = "MCI"
        elif title_folded.startswith("ad patient"):
            endpoint = "AD"
        elif title_folded.startswith("cognitively healthy control"):
            endpoint = "HC"
        else:
            raise ValueError(f"unrecognized official sample title category: {title}")
        rows.append(
            {
                "assessment_index": str(index),
                "endpoint": endpoint,
                "sex": sex,
                "age_reported": "no",
                "specimen": "blood_leukocytes",
            }
        )
    return rows


def audit() -> tuple[list[dict[str, str]], dict[str, object]]:
    rows = parse_metadata(HEADER_PATH)
    endpoint_counts = {endpoint: sum(row["endpoint"] == endpoint for row in rows) for endpoint in ("AD", "HC", "MCI")}
    sex_counts = {sex: sum(row["sex"] == sex for row in rows) for sex in sorted({row["sex"] for row in rows})}
    beta_available = BETA_PATH.exists() and BETA_PATH.stat().st_size == OFFICIAL_BETA_EXPECTED_BYTES
    metadata_status = (
        "METADATA_ONLY_CROSS_OMICS_AUDIT_SOURCE_LOCKED_SEPARATE_ASSOCIATION_AVAILABLE"
        if beta_available
        else "METADATA_ONLY_CROSS_OMICS_AUDIT_BETA_NOT_AVAILABLE"
    )
    beta_retrieval_status = (
        "SOURCE_LOCKED_SEPARATE_ASSOCIATION_AVAILABLE"
        if beta_available
        else "NOT_RETRIEVED_CONNECTION_FAILURE_OR_TIMEOUT"
    )
    result = {
        "status": metadata_status,
        "accession": "GSE208623",
        "official_series_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208623",
        "official_platform": "GPL21145",
        "assay": "Illumina Infinium MethylationEPIC DNA methylation array",
        "specimen": "peripheral blood leukocytes",
        "cohort_country": "China",
        "metadata_rows": len(rows),
        "endpoint_counts": endpoint_counts,
        "sex_available": all(row["sex"] in {"female", "male"} for row in rows),
        "sex_counts": sex_counts,
        "age_available": False,
        "age_metadata_status": "NOT_REPORTED_IN_INSPECTED_OFFICIAL_SERIES_MATRIX_METADATA",
        "metadata_sha256": sha256_file(HEADER_PATH),
        "processed_beta_source_url": BETA_SOURCE_URL,
        "processed_beta_expected_bytes": OFFICIAL_BETA_EXPECTED_BYTES,
        "beta_matrix_available": beta_available,
        "beta_matrix_sha256": sha256_file(BETA_PATH) if beta_available else None,
        "beta_retrieval_status": beta_retrieval_status,
        "metadata_module_performed_cpg_association": False,
        "metadata_module_performed_methylation_expression_overlap": False,
        "separate_source_locked_association_available": beta_available,
        "mRNA_model_fitting": False,
        "mRNA_model_tuning_or_calibration": False,
        "mRNA_model_validation": False,
        "mr_or_causal_analysis_performed": False,
        "claim_boundary": (
            "Independent blood-leukocyte DNA-methylation cohort metadata only. A separate source-locked "
            "sex-adjusted association audit is available; this metadata module did not perform it. Not mRNA "
            "model validation, model fitting, tuning, calibration, pooled analysis, methylation-expression "
            "mediation, cell-composition proof, or causal evidence."
        ),
    }
    return rows, result


def main(output_dir: Path = OUTPUT_DIR) -> int:
    rows, result = audit()
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "gse208623_metadata_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "accession",
            "metadata_rows",
            "endpoint_counts",
            "sex_counts",
            "sex_available",
            "age_available",
            "metadata_sha256",
        ]
        aggregate_row = {
            "accession": result["accession"],
            "metadata_rows": result["metadata_rows"],
            "endpoint_counts": json.dumps(result["endpoint_counts"], sort_keys=True),
            "sex_counts": json.dumps(result["sex_counts"], sort_keys=True),
            "sex_available": result["sex_available"],
            "age_available": result["age_available"],
            "metadata_sha256": result["metadata_sha256"],
        }
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(aggregate_row)
    (output_dir / "gse208623_metadata_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
