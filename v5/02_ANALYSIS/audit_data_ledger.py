from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from run_pairwise_external import (
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    GSE140829_MATRIX,
    GPL15988_FULL,
    load_dev,
    load_external,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "01_DATA_AUDIT"
GSE249477_RAW = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_raw_count_normalize_04-10-2025.csv.gz"
)
GSE249477_SERIES = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_series_matrix.txt.gz"
)
GSE282742_AUDIT = Path(
    r"E:\20241004_MCI\20260709\external_public_data_20260721_151749\GSE282742_sample_subject_audit.csv"
)
GSE282742_MATRIX = Path(
    r"E:\20241004_MCI\20260709\external_public_data_20260721_151749\GSE282742_Expected_count.full.gz"
)
GSE108215_AUDIT = ROOT / "03_RESULTS" / "gse108215_paired_audit.json"
GSE1408_AUDIT = ROOT / "03_RESULTS" / "gse1408_overlap_audit.json"
GSE108215_MATRIX = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse108215_task10\GSE108215_series_matrix.complete.txt.gz"
)
GSE108215_ANNOTATION = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse108215_task10\GPL570.annot.complete.gz"
)


def parse_gse249477() -> dict[str, object]:
    titles: list[str] = []
    with gzip.open(GSE249477_SERIES, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:]
                break
    counts = {"HC": 0, "MCI": 0, "AD": 0, "OTHER": 0}
    for title in titles:
        lower = title.casefold()
        if "mild cognitive impairment" in lower:
            counts["MCI"] += 1
        elif "cognitively normal" in lower:
            counts["HC"] += 1
        elif "alzheimer" in lower:
            counts["AD"] += 1
        else:
            counts["OTHER"] += 1
    return {
        "dataset": "GSE249477",
        "role": "cross-platform stress test only",
        "assay": "blood bulk RNA-seq",
        "specimen": "peripheral blood",
        "raw_n": len(titles),
        "eligible_n": counts["HC"] + counts["MCI"] + counts["AD"],
        "hc_n": counts["HC"],
        "mci_n": counts["MCI"],
        "ad_n": counts["AD"],
        "analysis_unit": "sample",
        "source_paths": f"{GSE249477_RAW};{GSE249477_SERIES}",
    }


def parse_gse282742() -> dict[str, object]:
    data = pd.read_csv(GSE282742_AUDIT)
    counts = data["group"].value_counts().to_dict()
    return {
        "dataset": "GSE282742",
        "role": "separate progression branch",
        "assay": "blood bulk RNA-seq",
        "specimen": "peripheral blood",
        "raw_n": int(len(data)),
        "eligible_n": int(len(data)),
        "hc_n": 0,
        "mci_n": int(counts.get("P-MCI", 0) + counts.get("S-MCI", 0)),
        "ad_n": int(counts.get("AD", 0)),
        "p_mci_n": int(counts.get("P-MCI", 0)),
        "s_mci_n": int(counts.get("S-MCI", 0)),
        "unique_subjects": int(data["subject_id"].nunique()),
        "baseline_unique_subject_n": 52,
        "p_mci_to_ad_pairs": 11,
        "analysis_unit": "sample for cohort audit; unique subject representative record for progression association (no verified visit date)",
        "source_paths": str(GSE282742_AUDIT),
    }


def parse_gse108215() -> dict[str, object]:
    audit = json.loads(GSE108215_AUDIT.read_text(encoding="utf-8"))
    if audit.get("status") != "EXPLORATORY_PAIRED_INTERVENTION_RESPONSE_ONLY":
        raise RuntimeError("GSE108215 paired audit is not in its required exploratory-only state")
    paired = audit.get("paired_subject_counts", {})
    mapping = audit.get("mapping", {})
    if paired != {"HEP": 18, "MAP": 17} or mapping.get("mapped") != 10:
        raise RuntimeError("GSE108215 paired count or mapping boundary changed")
    return {
        "dataset": "GSE108215",
        "role": "supplementary paired intervention-response audit only",
        "assay": "whole-blood bulk microarray",
        "specimen": "peripheral whole blood",
        "raw_n": int(audit["sample_counts"]["total_expression_rows"]),
        "eligible_n": int(sum(paired.values())),
        "hc_n": 0,
        "mci_n": int(audit["sample_counts"]["total_expression_rows"]),
        "ad_n": 0,
        "analysis_unit": "paired subject-level post-minus-pre change; not a stage-validation unit",
        "source_paths": f"{GSE108215_MATRIX};{GSE108215_ANNOTATION}",
    }


def main() -> int:
    if not GSE1408_AUDIT.exists():
        raise RuntimeError("GSE140829/830/831 overlap audit must run before rebuilding the ledger")
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    common_dev = sorted(set(dev60_x.columns).intersection(dev61_x.columns))
    ext_x, ext_meta, ext_audit = load_external(GSE140829_MATRIX, GPL15988_FULL, set(common_dev))
    rows = [
        {
            "dataset": "GSE63060",
            "role": "development; AddNeuroMed family",
            "assay": "blood bulk microarray",
            "specimen": "peripheral whole blood",
            "raw_n": int(len(dev60_meta)),
            "eligible_n": int(len(dev60_meta)),
            "hc_n": int((dev60_meta["label"] == "HC").sum()),
            "mci_n": int((dev60_meta["label"] == "MCI").sum()),
            "ad_n": int((dev60_meta["label"] == "AD").sum()),
            "analysis_unit": "sample",
            "source_paths": f"{GSE60_MATRIX};{GSE60_LABELS}",
        },
        {
            "dataset": "GSE63061",
            "role": "development; AddNeuroMed family",
            "assay": "blood bulk microarray",
            "specimen": "peripheral whole blood",
            "raw_n": int(len(dev61_meta)),
            "eligible_n": int(len(dev61_meta)),
            "hc_n": int((dev61_meta["label"] == "HC").sum()),
            "mci_n": int((dev61_meta["label"] == "MCI").sum()),
            "ad_n": int((dev61_meta["label"] == "AD").sum()),
            "analysis_unit": "sample",
            "source_paths": f"{GSE61_MATRIX};{GSE61_LABELS}",
        },
        {
            "dataset": "GSE140829",
            "role": "external transport assessment; subject map unresolved",
            "assay": "blood bulk microarray",
            "specimen": "peripheral whole blood",
            "raw_n": int(ext_audit["series_sample_count"]),
            "eligible_n": int(ext_audit["retained_sample_count"]),
            "hc_n": int(ext_meta["label"].eq("HC").sum()),
            "mci_n": int(ext_meta["label"].eq("MCI").sum()),
            "ad_n": int(ext_meta["label"].eq("AD").sum()),
            "analysis_unit": "sample-row audit; subject independence unresolved",
            "source_paths": f"{GSE140829_MATRIX};{GPL15988_FULL}",
        },
        parse_gse249477(),
        parse_gse282742(),
        parse_gse108215(),
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    ledger = pd.DataFrame(rows)
    ledger.to_csv(OUT / "data_ledger.csv", index=False)
    hashes = {}
    for path in [
        GSE60_MATRIX,
        GSE60_LABELS,
        GSE61_MATRIX,
        GSE61_LABELS,
        GSE140829_MATRIX,
        GPL15988_FULL,
        GSE249477_RAW,
        GSE249477_SERIES,
        GSE282742_AUDIT,
        GSE282742_MATRIX,
        GSE108215_AUDIT,
        GSE108215_MATRIX,
        GSE108215_ANNOTATION,
    ]:
        if path.exists():
            hashes[str(path)] = sha256_file(path)
    audit = {
        "status": "AUDITED",
        "ledger": "data_ledger.csv",
        "common_dev_exact_symbol_gene_count": len(common_dev),
        "gse140829_platform_audit": ext_audit,
        "gse249477_official_title_labels": True,
        "gse282742_progression_sample_and_subject_units_separated": True,
        "gse108215_paired_intervention_audit_only": True,
        "source_hashes": hashes,
    }
    (OUT / "data_ledger_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = [
        "# MCI public-data ledger",
        "",
        "This ledger separates sample-level cohort counts from subject-level progression counts.",
        "It is generated from official public matrices/metadata and is the source for figure and manuscript n values.",
        "",
        ledger.to_markdown(index=False),
        "",
        "## Locked corrections",
        "",
        "- GSE249477: official series titles yield HC=21, MCI=20, AD=21 (n=62).",
        "- GSE282742: 116 sample rows correspond to 86 unique subjects; the P-MCI/S-MCI branch uses one metadata-order representative record per subject (n=52; 22 P-MCI, 30 S-MCI), not a verified baseline visit. Age/sex-adjusted analysis has n=51 (22 P-MCI, 29 S-MCI); 11 P-MCI-to-AD subject pairs are a separate sensitivity analysis.",
        "- GSE108215: 86 expression rows include 35 traceable paired MCI subjects (HEP=18, MAP=17); it is supplementary intervention-response context, not MCI-versus-HC validation.",
        "- No count is promoted to an external validation claim without assay, label, mapping and analysis-unit compatibility.",
    ]
    (OUT / "data_ledger.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(ledger.to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
