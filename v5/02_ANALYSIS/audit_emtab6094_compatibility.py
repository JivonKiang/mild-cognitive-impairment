"""Record the official BioStudies compatibility audit for E-MTAB-6094."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SOURCE_OUT = ROOT / "04_FIGURES" / "source_data"


def main() -> int:
    audit = {
        "status": "EXCLUDED_FROM_PRIMARY_MCI_HC",
        "accession": "E-MTAB-6094",
        "title": "Whole genome transcriptional analysis of Alzheimer's disease patients compared with age-matched non-demented controls",
        "assay": "transcription profiling by array; Agilent array",
        "organism": "Homo sapiens",
        "organism_part": "blood",
        "cell_type": "mononuclear cell",
        "specimen_interpretation": "PBMC; not peripheral whole blood",
        "biostudies_factor_sample_count": 35,
        "factor_counts": {"normal": 13, "AD_CDR1": 15, "AD_CDR3": 7},
        "mci_factor_count": 0,
        "description_claimed_ad_count": 25,
        "description_claimed_control_count": 15,
        "metadata_count_conflict": True,
        "subject_level_map_verified": False,
        "grouped_expression_analysis_performed": False,
        "decision": "Exclude from primary MCI-versus-HC endpoint: no MCI group, PBMC rather than whole blood, and a study-description versus factor-table count conflict.",
        "allowed_role": "orthogonal PBMC AD-versus-control context after SDRF reconciliation",
        "forbidden_claims": ["independent MCI-versus-HC validation", "whole-blood replication", "clinical diagnostic validation"],
        "official_source": "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-6094",
        "record_page": "https://www.ebi.ac.uk/biostudies/studies/E-MTAB-6094",
        "access_date": "2026-08-13",
        "source_response_persisted": False,
        "source_response_note": "The official API response was inspected on the access date. A subsequent raw-response save was interrupted by transient network refusal; no raw response hash is claimed.",
    }
    payload = json.dumps(audit, sort_keys=True, ensure_ascii=True).encode("utf-8")
    audit["audit_payload_sha256"] = hashlib.sha256(payload).hexdigest()
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    for destination in [OUT / "emtab6094_compatibility_audit.json", SOURCE_OUT / "emtab6094_compatibility_audit.json"]:
        destination.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
