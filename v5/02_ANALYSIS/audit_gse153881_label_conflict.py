"""Audit the GSE153881 public label conflict without performing group analysis."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_BUNDLE = Path(r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse153881_task10")
MATRIX = RUN_BUNDLE / "GSE153881_series_matrix.complete.txt.gz"
ARCHIVE = RUN_BUNDLE / "GSE153881_RAW.tar"
OUT = ROOT / "03_RESULTS"
SOURCE_OUT = ROOT / "04_FIGURES" / "source_data"
FILENAME_PATTERN = re.compile(r"^(GSM\d+)_.*_(healthy|MCI|AD)_.*_geneCOUNT\.txt\.gz$", re.IGNORECASE)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_label(value: str) -> str:
    lower = value.lower()
    if "healthy" in lower:
        return "HC"
    if "mild cognitive impairment" in lower or lower == "mci":
        return "MCI"
    if "alzheimer" in lower or lower == "ad":
        return "AD"
    raise ValueError(f"unrecognized diagnosis label: {value}")


def metadata_labels() -> dict[str, str]:
    accessions: list[str] = []
    diagnoses: list[str] = []
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_characteristics_ch1") and "diagnosis:" in line.lower():
                diagnoses = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:]
    if len(accessions) != 15 or len(diagnoses) != 15:
        raise RuntimeError("GSE153881 metadata lacks 15 accessions and 15 diagnosis records")
    result = {}
    for accession, diagnosis in zip(accessions, diagnoses):
        result[accession] = normalize_label(diagnosis.split(":", 1)[-1].strip())
    return result


def archive_labels() -> tuple[dict[str, str], int, int]:
    labels: dict[str, str] = {}
    valid_count_files = 0
    parseable_count_files = 0
    with tarfile.open(ARCHIVE, "r") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        for member in members:
            match = FILENAME_PATTERN.fullmatch(Path(member.name).name)
            if match is None:
                raise RuntimeError(f"unexpected GSE153881 archive member: {member.name}")
            accession, label = match.groups()
            labels[accession] = normalize_label(label)
            valid_count_files += 1
            handle = archive.extractfile(member)
            if handle is None:
                raise RuntimeError(f"unable to read archive member: {member.name}")
            with gzip.open(handle, "rt", encoding="utf-8", errors="replace") as count_file:
                header = count_file.readline().rstrip("\r\n").split("\t")
                first_data = count_file.readline().rstrip("\r\n").split("\t")
            if len(header) == 2 and len(first_data) == 2 and header[0] == "Tracking_ID":
                parseable_count_files += 1
    return labels, valid_count_files, parseable_count_files


def main() -> int:
    for path in [MATRIX, ARCHIVE]:
        if not path.exists():
            raise FileNotFoundError(f"required source missing: {path}")
    geo = metadata_labels()
    raw, count_files, parseable_count_files = archive_labels()
    if set(geo) != set(raw) or count_files != 15 or parseable_count_files != 15:
        raise RuntimeError("GSE153881 accession coverage or count-file parse gate failed")
    conflicts = [
        {"accession": accession, "geo_matrix_label": geo[accession], "raw_filename_label": raw[accession]}
        for accession in sorted(geo)
        if geo[accession] != raw[accession]
    ]
    audit = {
        "status": "BLOCKED_PENDING_LABEL_RECONCILIATION",
        "dataset": "GSE153881",
        "assay": "plasma extracellular-vesicle RNA-seq",
        "specimen": "plasma extracellular vesicles",
        "supplementary_count_files": count_files,
        "parseable_count_files": parseable_count_files,
        "metadata_sample_count": len(geo),
        "accession_label_conflict_count": len(conflicts),
        "conflict_summary": "GEO series-matrix diagnosis records and supplementary geneCOUNT filename labels disagree for MCI and AD accessions.",
        "grouped_expression_analysis_performed": False,
        "allowed_claim": "The public source records contain an accession-level label conflict requiring independent reconciliation before groupwise analysis.",
        "forbidden_claim": "No MCI-versus-HC, MCI-versus-AD, diagnostic, prognostic, mechanistic, or treatment result may be reported from this dataset at this stage.",
        "source_hashes": {str(MATRIX): sha256(MATRIX), str(ARCHIVE): sha256(ARCHIVE)},
        "official_sources": {
            "geo_record": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153881",
            "linked_publication": "https://pubmed.ncbi.nlm.nih.gov/33304899/",
            "access_date": "2026-08-13",
        },
        "conflicts": conflicts,
    }
    if len(conflicts) != 10:
        raise RuntimeError(f"expected 10 GSE153881 MCI/AD label conflicts, observed {len(conflicts)}")
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    for destination in [OUT / "gse153881_label_conflict_audit.json", SOURCE_OUT / "gse153881_label_conflict_audit.json"]:
        destination.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
