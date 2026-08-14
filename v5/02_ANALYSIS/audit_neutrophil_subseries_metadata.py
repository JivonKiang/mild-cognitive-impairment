"""Audit only the official metadata boundaries of APOE neutrophil SubSeries.

This module deliberately does not load expression values, fit models, or report
MCI-HC effects. It records why the SubSeries cannot validate a whole-blood model.
"""
from __future__ import annotations

import csv
import argparse
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "01_DATA_AUDIT" / "neutrophil_subseries"
OUTPUT_DIR = ROOT / "03_RESULTS"
ACCESSIONS = ("GSE243743", "GSE243744", "GSE243745")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_matrix_metadata(path: Path) -> list[dict[str, str]]:
    fields: dict[str, list[str]] = {}
    characteristic_rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                if line.startswith("!series_matrix_table_begin"):
                    break
                continue
            key, *values = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            if key == "!Sample_characteristics_ch1":
                characteristic_rows.append(values)
            else:
                fields[key] = values
    gsm = fields.get("!Sample_geo_accession", [])
    titles = fields.get("!Sample_title", [])
    if not gsm or len(gsm) != len(titles) or any(len(values) != len(gsm) for values in characteristic_rows):
        raise ValueError(f"invalid or incomplete official sample metadata: {path.name}")
    rows = [{"gsm": gsm[index], "title": titles[index]} for index in range(len(gsm))]
    for values in characteristic_rows:
        field_name = values[0].split(":", 1)[0].strip().casefold().replace(" ", "_")
        if not field_name:
            raise ValueError(f"unparseable characteristic field in {path.name}")
        for index, value in enumerate(values):
            prefix, separator, content = value.partition(":")
            if not separator:
                raise ValueError(f"unparseable characteristic value in {path.name}")
            rows[index][field_name] = content.strip()
    return rows


def audit() -> tuple[list[dict[str, object]], dict[str, object]]:
    all_gsm: dict[str, set[str]] = {}
    rows: list[dict[str, object]] = []
    for accession in ACCESSIONS:
        path = INPUT_DIR / f"{accession}_series_matrix.txt.gz"
        sample_rows = parse_matrix_metadata(path)
        all_gsm[accession] = {str(row["gsm"]) for row in sample_rows}
        disease_counts: dict[str, int] = {}
        for row in sample_rows:
            disease = str(row.get("disease_state", "MISSING"))
            disease_counts[disease] = disease_counts.get(disease, 0) + 1
        rows.append(
            {
                "accession": accession,
                "official_metadata_rows": len(sample_rows),
                "gsm_unique": len(all_gsm[accession]),
                "disease_counts": json.dumps(disease_counts, sort_keys=True),
                "sex_values": json.dumps(sorted({str(row.get("sex", "MISSING")) for row in sample_rows})),
                "genotype_values": json.dumps(sorted({str(row.get("genotype", "MISSING")) for row in sample_rows})),
                "selection_values": json.dumps(sorted({str(row.get("il18r1_profile", "NOT_REPORTED")) for row in sample_rows})),
                "metadata_sha256": sha256_file(path),
                "allowed_role": "directional purified-neutrophil context only",
                "forbidden_role": "whole-blood validation, pooled cohort, independent replication, causal mechanism",
            }
        )
    pairwise_overlap = {
        f"{left}_{right}": len(all_gsm[left] & all_gsm[right])
        for index, left in enumerate(ACCESSIONS)
        for right in ACCESSIONS[index + 1 :]
    }
    result = {
        "status": "METADATA_ONLY_NEUTROPHIL_CONTEXT_AUDIT",
        "accessions": list(ACCESSIONS),
        "pairwise_gsm_overlap": pairwise_overlap,
        "subject_identifier_status": "NOT_AVAILABLE_IN_OFFICIAL_SERIES_MATRIX_METADATA",
        "subject_independence_established": False,
        "expression_loaded": False,
        "predictive_model_fitting": False,
        "claim_boundary": "Purified-neutrophil, APOE/sex/selection-restricted metadata context only; not whole-blood validation, independent replication, pooled analysis, or causal mechanism.",
    }
    return rows, result


def main(output_dir: Path = OUTPUT_DIR) -> int:
    rows, result = audit()
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "neutrophil_subseries_metadata_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output_dir / "neutrophil_subseries_metadata_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
