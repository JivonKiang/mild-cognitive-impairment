"""Audit the count-column to GEO-sample join for GSE270454.

The frozen score audit previously assumed that the first 45 GEO metadata
records were in the same order as the count matrix columns.  This audit
rebuilds the join using the group and ordinal encoded in the count-column
names and the official GEO sample metadata.  It writes a local, traceable
mapping plus a sanitized aggregate audit for manuscript-facing provenance.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


PROJECT = Path(r"E:\20260805 Biomni\MCI_PUBLIC_RELIABILITY_20260812_V5")
DATA_ROOT = Path(r"E:\20241004_MCI\20260805_Human_Genomics\03_REPRODUCIBILITY")
COUNTS = DATA_ROOT / "GSE270454_RNAseq-combined-counts-matrix.csv.gz"
SERIES = DATA_ROOT / "GSE270454_series_matrix.txt.gz"
LOCAL_MAPPING = PROJECT / "01_DATA_AUDIT" / "gse270454_join_mapping_local.csv"
AUDIT_JSON = PROJECT / "03_RESULTS" / "gse270454_join_audit.json"
SUMMARY_CSV = PROJECT / "04_FIGURES" / "source_data" / "gse270454_join_summary.csv"

COUNT_PATTERN = re.compile(r"^(ASM|ASO|MCI|AD)_(\d+)$", re.I)
TOKEN_PATTERN = re.compile(r"\b(ASM|ASO|MCI|AD)[ _-]?(\d+)\b", re.I)
GROUP_PATTERN = re.compile(r"\b(ASM|ASO|MCI|AD)\b", re.I)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gse_sample_rows(path: Path) -> list[dict[str, object]]:
    """Parse official GEO sample rows, preserving repeated characteristics."""
    columns: dict[str, list[str]] = {}
    repeated_characteristics: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw in handle:
            if not raw.startswith("!Sample_"):
                continue
            fields = next(csv.reader([raw.rstrip("\n")], delimiter="\t"), [])
            if not fields:
                continue
            key = fields[0]
            values = fields[1:]
            if key == "!Sample_characteristics_ch1":
                repeated_characteristics.append(values)
            elif key not in columns:
                columns[key] = values

    accessions = columns.get("!Sample_geo_accession", [])
    titles = columns.get("!Sample_title", [])
    sources = columns.get("!Sample_source_name_ch1", [])
    n = len(accessions)
    records: list[dict[str, object]] = []
    for index in range(n):
        characteristics: list[str] = []
        for row in repeated_characteristics:
            if index < len(row) and row[index].strip():
                characteristics.append(row[index].strip())
        records.append(
            {
                "geo_index": index + 1,
                "accession": accessions[index].strip() if index < len(accessions) else "",
                "title": titles[index].strip() if index < len(titles) else "",
                "source": sources[index].strip() if index < len(sources) else "",
                "characteristics": characteristics,
            }
        )
    return records


def parse_characteristics(values: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        if ":" not in value:
            continue
        key, item = value.split(":", 1)
        parsed[key.strip().lower()] = item.strip()
    return parsed


def normalize_group(value: str) -> str | None:
    matches = {match.group(1).upper() for match in GROUP_PATTERN.finditer(value)}
    if len(matches) == 1:
        return next(iter(matches))
    return None


def extract_group_ordinal(record: dict[str, object]) -> tuple[str | None, int | None, list[str]]:
    title = str(record.get("title", "")).strip()
    title_lower = title.lower()
    semantic_group: str | None = None
    # GEO titles use descriptive cohort labels rather than the abbreviations
    # used by the supplementary count matrix.  The longer label must be
    # checked first because it contains the phrase "older adult".
    if "middle-aged/older adult" in title_lower:
        semantic_group = "ASM"
    elif "older adult" in title_lower:
        semantic_group = "ASO"
    elif "mci patient" in title_lower:
        semantic_group = "MCI"
    elif "ad patient" in title_lower:
        semantic_group = "AD"

    values = [
        title,
        str(record.get("source", "")),
        *[str(item) for item in record.get("characteristics", [])],
    ]
    text = " | ".join(values)
    tokens = [(match.group(1).upper(), int(match.group(2))) for match in TOKEN_PATTERN.finditer(text)]
    unique_tokens = list(dict.fromkeys(tokens))
    groups = sorted({group for group, _ in unique_tokens})
    ordinals = sorted({ordinal for _, ordinal in unique_tokens})
    conflicts: list[str] = []
    if semantic_group is None:
        if len(groups) != 1:
            conflicts.append("group_token_conflict" if groups else "group_token_missing")
        group = groups[0] if len(groups) == 1 else None
    else:
        group = semantic_group
        if groups and groups != [semantic_group]:
            conflicts.append("group_semantic_token_conflict")

    # For descriptive GEO titles, the terminal integer is the sample ordinal;
    # do not mistake age values in characteristics for sample identifiers.
    ordinal_match = re.search(r",\s*(\d+)\s*$", title)
    if ordinal_match:
        ordinal = int(ordinal_match.group(1))
        if ordinals and ordinals != [ordinal]:
            conflicts.append("ordinal_title_token_conflict")
    elif len(ordinals) == 1:
        ordinal = ordinals[0]
    else:
        ordinal = None
        conflicts.append("ordinal_token_conflict" if ordinals else "ordinal_token_missing")
    return group, ordinal, conflicts


def parse_count_columns(path: Path) -> list[dict[str, object]]:
    header = pd.read_csv(path, compression="gzip", nrows=0).columns.tolist()
    records: list[dict[str, object]] = []
    for column in header[1:]:
        name = str(column).strip()
        match = COUNT_PATTERN.fullmatch(name)
        if match:
            records.append(
                {
                    "count_column": name,
                    "group": match.group(1).upper(),
                    "ordinal": int(match.group(2)),
                }
            )
    return records


def build_mapping(count_records: list[dict[str, object]], geo_records: list[dict[str, object]]) -> tuple[pd.DataFrame, dict[str, object]]:
    enriched: list[dict[str, object]] = []
    for record in geo_records:
        characteristics = record["characteristics"]
        assert isinstance(characteristics, list)
        parsed = parse_characteristics(characteristics)
        group, ordinal, conflicts = extract_group_ordinal(record)
        tissue_text = " | ".join(
            [str(record.get("source", "")), *[str(item) for item in characteristics]]
        ).lower()
        enriched.append(
            {
                **record,
                "group": group,
                "ordinal": ordinal,
                "conflicts": conflicts,
                "age": parsed.get("age", ""),
                "sex": parsed.get("sex", ""),
                "tissue_text": tissue_text,
                "whole_blood": "whole blood" in tissue_text,
            }
        )

    rna_records = [record for record in enriched if record["whole_blood"]]
    lookup: dict[tuple[str, int], list[dict[str, object]]] = {}
    for record in rna_records:
        group = record["group"]
        ordinal = record["ordinal"]
        if group is not None and ordinal is not None:
            lookup.setdefault((str(group), int(ordinal)), []).append(record)

    rows: list[dict[str, object]] = []
    for count_record in count_records:
        key = (str(count_record["group"]), int(count_record["ordinal"]))
        candidates = lookup.get(key, [])
        match = candidates[0] if len(candidates) == 1 else None
        rows.append(
            {
                "count_column": count_record["count_column"],
                "count_group": key[0],
                "count_ordinal": key[1],
                "candidate_count": len(candidates),
                "geo_accession": match["accession"] if match else "",
                "geo_index": match["geo_index"] if match else "",
                "title": match["title"] if match else "",
                "source": match["source"] if match else "",
                "geo_group": match["group"] if match else "",
                "geo_ordinal": match["ordinal"] if match else "",
                "age": match["age"] if match else "",
                "sex": match["sex"] if match else "",
                "whole_blood": match["whole_blood"] if match else False,
                "metadata_conflicts": ";".join(match["conflicts"]) if match else "unresolved",
                "join_status": "unique" if match else ("ambiguous" if candidates else "unmatched"),
            }
        )
    mapping = pd.DataFrame(rows)
    counts_by_group = {
        group: int(sum(str(row["count_group"]) == group for row in rows))
        for group in ("ASM", "ASO", "MCI", "AD")
    }
    matched_by_group = {
        group: int(
            sum(
                str(row["count_group"]) == group and row["join_status"] == "unique"
                for row in rows
            )
        )
        for group in ("ASM", "ASO", "MCI", "AD")
    }
    audit = {
        "status": "PASS" if len(mapping) == 45 and (mapping["join_status"] == "unique").all() else "FAIL",
        "dataset": "GSE270454",
        "join_method": "count-column group and ordinal to unique GEO title/source/characteristics group and ordinal token",
        "count_columns": len(count_records),
        "geo_sample_records": len(geo_records),
        "geo_whole_blood_records": len(rna_records),
        "unique_join_count": int((mapping["join_status"] == "unique").sum()),
        "unmatched_count": int((mapping["join_status"] == "unmatched").sum()),
        "ambiguous_count": int((mapping["join_status"] == "ambiguous").sum()),
        "count_columns_by_group": counts_by_group,
        "unique_joins_by_group": matched_by_group,
        "geo_group_token_counts": {
            group: int(sum(record["group"] == group for record in rna_records))
            for group in ("ASM", "ASO", "MCI", "AD")
        },
        "group_or_ordinal_conflict_records": int(sum(bool(record["conflicts"]) for record in rna_records)),
        "all_matched_source_whole_blood": bool(
            len(mapping) > 0 and mapping.loc[mapping["join_status"] == "unique", "whole_blood"].all()
        ),
        "whole_blood_true_count": int(mapping["whole_blood"].sum()),
        "age_nonempty_unique_joins": int(
            mapping.loc[mapping["join_status"] == "unique", "age"].astype(str).str.strip().ne("").sum()
        ),
        "sex_nonempty_unique_joins": int(
            mapping.loc[mapping["join_status"] == "unique", "sex"].astype(str).str.strip().ne("").sum()
        ),
        "source_hashes": {
            "counts": sha256(COUNTS),
            "series_matrix": sha256(SERIES),
        },
        "sample_level_identifiers_in_formal_summary": False,
    }
    return mapping, audit


def main() -> None:
    count_records = parse_count_columns(COUNTS)
    geo_records = parse_gse_sample_rows(SERIES)
    mapping, audit = build_mapping(count_records, geo_records)

    LOCAL_MAPPING.parent.mkdir(parents=True, exist_ok=True)
    AUDIT_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(LOCAL_MAPPING, index=False)
    AUDIT_JSON.write_text(json.dumps(audit, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")

    summary = pd.DataFrame(
        [
            {"metric": "status", "value": audit["status"]},
            {"metric": "count_columns", "value": audit["count_columns"]},
            {"metric": "geo_sample_records", "value": audit["geo_sample_records"]},
            {"metric": "unique_join_count", "value": audit["unique_join_count"]},
            {"metric": "unmatched_count", "value": audit["unmatched_count"]},
            {"metric": "ambiguous_count", "value": audit["ambiguous_count"]},
            {"metric": "whole_blood_true_count", "value": audit["whole_blood_true_count"]},
            {"metric": "age_nonempty_unique_joins", "value": audit["age_nonempty_unique_joins"]},
            {"metric": "sex_nonempty_unique_joins", "value": audit["sex_nonempty_unique_joins"]},
        ]
    )
    summary.to_csv(SUMMARY_CSV, index=False)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    if audit["status"] != "PASS":
        raise SystemExit("GSE270454 count-to-GEO join audit failed")


if __name__ == "__main__":
    main()
