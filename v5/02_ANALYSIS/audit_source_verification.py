from __future__ import annotations

"""Materialize a conservative, offline-verifiable cohort-paper provenance audit."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "03_RESULTS"
ACCESS_DATE = "2026-08-13"
GSE140829_PREPRINT_DOI = "10.1101/2019.12.13.875112"
GSE140829_PREPRINT_TITLE = (
    "Systems-level analysis of peripheral blood gene expression in dementia patients "
    "reveals an innate immune response shared across multiple disorders"
)
GSE140829_CROSSREF_URL = f"https://api.crossref.org/works/{GSE140829_PREPRINT_DOI}"


RECORDS = [
    {
        "accession": "GSE63060",
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060",
        "local_series_metadata": r"E:\20241004_MCI\20260709\data\GSE63060_series_matrix.txt.gz",
        "geo_linked_pmid": "26343147",
        "expected_title": "A novel multi-tissue RNA diagnostic of healthy ageing relates to cognitive health status.",
        "dataset_role_boundary": "AddNeuroMed development-family batch 1; not independent external validation.",
        "matching_status": "ACCESSION_AND_LINKED_PAPER_VERIFIED",
    },
    {
        "accession": "GSE63061",
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061",
        "local_series_metadata": r"E:\20241004_MCI\20260709\data\GSE63061_series_matrix.txt.gz",
        "geo_linked_pmid": "26343147",
        "expected_title": "A novel multi-tissue RNA diagnostic of healthy ageing relates to cognitive health status.",
        "dataset_role_boundary": "AddNeuroMed development-family batch 2; not independent external validation.",
        "matching_status": "ACCESSION_AND_LINKED_PAPER_VERIFIED",
    },
    {
        "accession": "GSE140829",
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829",
        "local_series_metadata": r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse140829_task8\GSE140829_series_matrix.fresh.txt.gz",
        "geo_linked_pmid": "",
        "expected_title": GSE140829_PREPRINT_TITLE,
        "preprint_doi": GSE140829_PREPRINT_DOI,
        "dataset_role_boundary": "Sample-row transport audit only; subject independence unresolved; excluded from pooled inference.",
        "matching_status": "ACCESSION_AND_TITLE_MATCHED_PREPRINT_VERIFIED",
    },
    {
        "accession": "GSE249477",
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249477",
        "local_series_metadata": r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_series_matrix.txt.gz",
        "geo_linked_pmid": "40034360",
        "expected_title": "Blood RNA transcripts show changes in inflammation and lipid metabolism in Alzheimer's disease and mitochondrial function in mild cognitive impairment.",
        "dataset_role_boundary": "Small independent RNA-seq cross-platform stress test; not positive replication or clinical validation.",
        "matching_status": "ACCESSION_AND_LINKED_PAPER_VERIFIED",
    },
    {
        "accession": "GSE108215",
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108215",
        "local_series_metadata": r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse108215_task10\GSE108215_series_matrix.complete.txt.gz",
        "geo_linked_pmid": "29900257",
        "expected_title": "Dataset on gene expression in the elderly after Mindfulness Awareness Practice or Health Education Program.",
        "dataset_role_boundary": "Supplementary paired intervention audit only; no HC comparator and not a primary endpoint validation cohort.",
        "matching_status": "ACCESSION_AND_LINKED_PAPER_VERIFIED",
    },
    {
        "accession": "GSE282742",
        "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282742",
        "local_series_metadata": "",
        "geo_linked_pmid": "39754267",
        "expected_title": "Transcriptomic predictors of rapid progression from mild cognitive impairment to Alzheimer's disease.",
        "dataset_role_boundary": "Separate progression-association branch; not independent MCI-vs-HC predictive validation and not causal evidence.",
        "matching_status": "LINKED_PAPER_VERIFIED_ACCESSION_HEADER_NOT_LOCAL",
    },
]
GSE140829_GDS_SUMMARY = RESULTS / "source_verification_gds_gse140829_esummary.json"
GSE140829_GDS_PUBMED_LINK = RESULTS / "source_verification_gse140829_gds_pubmed_elink.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pubmed_record(pmid: str) -> tuple[dict[str, object], str]:
    path = RESULTS / f"source_verification_pubmed_{pmid}.json"
    if not path.exists():
        raise FileNotFoundError(f"missing saved official PubMed response: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    return data["result"][pmid], sha256(path)


def doi(record: dict[str, object]) -> str:
    for item in record.get("articleids", []):
        if item.get("idtype") == "doi":
            return str(item["value"])
    return ""


def build_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for declared in RECORDS:
        row = dict(declared)
        pmid = row["geo_linked_pmid"]
        if pmid:
            record, response_hash = pubmed_record(pmid)
            title = str(record["title"])
            if title != row["expected_title"]:
                raise ValueError(f"PubMed title mismatch for {row['accession']}: {title!r}")
            row.update(
                {
                    "pubmed_title": title,
                    "pubmed_journal": str(record["fulljournalname"]),
                    "pubmed_pubdate": str(record["pubdate"]),
                    "pubmed_doi": doi(record),
                    "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "saved_pubmed_response_sha256": response_hash,
                    "publication_type": "journal_article",
                }
            )
        elif row.get("preprint_doi"):
            row.update(
                {
                    "pubmed_title": GSE140829_PREPRINT_TITLE,
                    "pubmed_journal": "bioRxiv",
                    "pubmed_pubdate": "2019-12-13",
                    "pubmed_doi": GSE140829_PREPRINT_DOI,
                    "pubmed_url": f"https://doi.org/{GSE140829_PREPRINT_DOI}",
                    "saved_pubmed_response_sha256": "",
                    "publication_type": "preprint",
                }
            )
        else:
            row.update({
                "pubmed_title": "",
                "pubmed_journal": "",
                "pubmed_pubdate": "",
                "pubmed_doi": "",
                "pubmed_url": "",
                "saved_pubmed_response_sha256": "",
                "publication_type": "",
            })
        local_path = Path(row["local_series_metadata"]) if row["local_series_metadata"] else None
        row["local_series_metadata_sha256"] = sha256(local_path) if local_path and local_path.exists() else ""
        row["local_series_metadata_status"] = "PRESENT" if local_path and local_path.exists() else "NOT_AVAILABLE_IN_LOCAL_AUDIT_INPUTS"
        row["access_date"] = ACCESS_DATE
        rows.append(row)
    return rows


def verify_gse140829_no_direct_pubmed_link() -> dict[str, object]:
    """Verify the saved official Entrez evidence for the deliberately pending link."""
    for path in [GSE140829_GDS_SUMMARY, GSE140829_GDS_PUBMED_LINK]:
        if not path.exists():
            raise FileNotFoundError(f"missing saved official GSE140829 response: {path}")
    summary = json.loads(GSE140829_GDS_SUMMARY.read_text(encoding="utf-8"))
    record = summary["result"]["200140829"]
    expected = {
        "accession": "GSE140829",
        "taxon": "Homo sapiens",
        "gpl": "15988",
        "gdstype": "Expression profiling by array",
    }
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(f"GSE140829 official GDS field changed: {key}")
    links = json.loads(GSE140829_GDS_PUBMED_LINK.read_text(encoding="utf-8"))
    linkset = links.get("linksets", [{}])[0]
    if linkset != {"dbfrom": "gds", "ids": ["200140829"]}:
        raise ValueError("GSE140829 official GDS-to-PubMed link response is no longer empty")
    return {
        "gds_summary_sha256": sha256(GSE140829_GDS_SUMMARY),
        "gds_to_pubmed_elink_sha256": sha256(GSE140829_GDS_PUBMED_LINK),
        "gds_to_pubmed_linked_pmids": [],
        "official_metadata": expected,
    }


def verify_gse140829_title_matched_preprint() -> dict[str, object]:
    """Verify an exact GEO-title to registered-preprint correspondence.

    This traceability link does not change the subject-independence, validation,
    biological, or peer-review evidence boundaries of the transport audit.
    """
    official = verify_gse140829_no_direct_pubmed_link()["official_metadata"]
    summary = json.loads(GSE140829_GDS_SUMMARY.read_text(encoding="utf-8"))
    geo_title = str(summary["result"]["200140829"]["title"]).removesuffix(" [ad_mci]")
    if geo_title != GSE140829_PREPRINT_TITLE:
        raise ValueError("GSE140829 GEO title does not exactly match the registered preprint title")
    return {
        "doi": GSE140829_PREPRINT_DOI,
        "doi_url": f"https://doi.org/{GSE140829_PREPRINT_DOI}",
        "crossref_url": GSE140829_CROSSREF_URL,
        "crossref_access_date": "2026-08-14",
        "publication_type": "preprint",
        "preprint_date": "2019-12-13",
        "title_exact_match": True,
        "geo_title": geo_title,
        "boundary": "Exact title-matched preprint provenance only; no peer-reviewed-paper claim, subject-independence resolution, or validation upgrade.",
        "official_metadata": official,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    gse140829_link_evidence = verify_gse140829_no_direct_pubmed_link()
    gse140829_preprint = verify_gse140829_title_matched_preprint()
    columns = [
        "accession", "geo_url", "geo_linked_pmid", "pubmed_title", "pubmed_journal",
        "pubmed_pubdate", "pubmed_doi", "pubmed_url", "matching_status", "publication_type",
        "dataset_role_boundary", "local_series_metadata", "local_series_metadata_status",
        "local_series_metadata_sha256", "saved_pubmed_response_sha256", "access_date",
    ]
    frame = pd.DataFrame(rows)[columns].sort_values("accession", kind="stable")
    csv_path = args.output_dir / "source_verification_audit.csv"
    json_path = args.output_dir / "source_verification_audit.json"
    frame.to_csv(csv_path, index=False, lineterminator="\n")
    payload = {
        "audit_type": "offline_reproducible_cohort_paper_provenance",
        "access_date": ACCESS_DATE,
        "record_count": len(frame),
        "verified_accession_and_paper_count": int(frame["matching_status"].str.startswith("ACCESSION_AND").sum()),
        "pending_paper_link_accessions": frame.loc[
            frame["matching_status"].eq("ACCESSION_METADATA_AUDITED_PAPER_LINK_PENDING"), "accession"
        ].tolist(),
        "gse140829_official_link_evidence": gse140829_link_evidence,
        "gse140829_title_matched_preprint": gse140829_preprint,
        "boundary": "This audit verifies metadata linkage and traceability only. It does not establish biological mechanism, causal effect, clinical utility, subject independence, or external predictive validation.",
        "csv_sha256": sha256(csv_path),
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
