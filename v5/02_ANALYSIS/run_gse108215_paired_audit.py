"""Supplementary paired audit for the public GSE108215 MCI intervention study.

This analysis is intentionally separate from MCI-versus-HC model development
and external transport assessment. It maps the frozen panel by exact GPL570
symbols, builds a partial score without arm/timepoint labels, and summarizes
only within-subject pre/post changes. No treatment outcome or diagnostic claim
is supported by this small, incomplete-follow-up cohort.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
RUN_BUNDLE = Path(r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse108215_task10")
MATRIX = RUN_BUNDLE / "GSE108215_series_matrix.complete.txt.gz"
ANNOTATION = RUN_BUNDLE / "GPL570.annot.complete.gz"
PANEL = Path(r"E:\20241004_MCI\20260709\mci_development_optimized_panel_20260722_161919\development_locked_12_gene_panel.csv")
OUT = ROOT / "03_RESULTS"
SOURCE_OUT = ROOT / "04_FIGURES" / "source_data"

TITLE_PATTERN = re.compile(
    r"blood,(?P<arm>HEP|MAP),(?P<time>Pre|Post)-treatment [EA]#\d+-(?P<subject>M\d+)-G"
)
EXPECTED_MISSING = ["PHF15", "LAT1-3TM"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_title(title: str) -> tuple[str, str, str]:
    """Return randomized arm, collection timepoint and local subject token."""
    match = TITLE_PATTERN.fullmatch(title.strip())
    if match is None:
        raise ValueError(f"unparseable GSE108215 title: {title}")
    return match.group("arm"), match.group("time"), match.group("subject")


def quoted_fields(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def read_series_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read GEO sample metadata and processed probe signal values."""
    metadata_fields: dict[str, list[str]] = {}
    table_marker_index: int | None = None
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            if line.startswith("!series_matrix_table_begin"):
                table_marker_index = index
                break
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
                fields = quoted_fields(line.rstrip("\n"))
                metadata_fields[fields[0]] = fields[1:]
    if table_marker_index is None:
        raise RuntimeError("GSE108215 matrix has no expression table")

    titles = metadata_fields.get("!Sample_title", [])
    accessions = metadata_fields.get("!Sample_geo_accession", [])
    if len(titles) != 86 or len(accessions) != 86 or len(titles) != len(accessions):
        raise RuntimeError("GSE108215 metadata does not contain the expected 86 matched sample rows")
    metadata = pd.DataFrame({"geo_accession": accessions, "title": titles})

    matrix = pd.read_csv(MATRIX, sep="\t", compression="gzip", skiprows=table_marker_index + 1)
    matrix.columns = [str(column).strip().strip('"') for column in matrix.columns]
    if "ID_REF" not in matrix.columns:
        raise RuntimeError("GSE108215 matrix lacks the ID_REF probe column")
    sample_columns = [column for column in matrix.columns if column != "ID_REF"]
    if sample_columns != accessions:
        raise RuntimeError("GSE108215 expression columns do not exactly match GEO metadata order")
    return metadata, matrix


def read_gpl570_annotation() -> pd.DataFrame:
    with gzip.open(ANNOTATION, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    try:
        start = next(index for index, line in enumerate(lines) if line.startswith("!platform_table_begin")) + 1
        end = next(index for index, line in enumerate(lines) if line.startswith("!platform_table_end"))
    except StopIteration as error:
        raise RuntimeError("GPL570 annotation table boundaries are absent") from error
    annotation = pd.read_csv(StringIO("".join(lines[start:end])), sep="\t", low_memory=False)
    if not {"ID", "Gene symbol"}.issubset(annotation.columns):
        raise RuntimeError("GPL570 annotation lacks ID and Gene symbol fields")
    return annotation[["ID", "Gene symbol"]].copy()


def exact_probe_map(annotation: pd.DataFrame, panel_genes: list[str]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for gene in panel_genes:
        matches = annotation["Gene symbol"].fillna("").astype(str).map(
            lambda value: gene in {item.strip().upper() for item in value.split("///")}
        )
        mapping[gene] = annotation.loc[matches, "ID"].astype(str).tolist()
    return mapping


def bootstrap_mean_ci(values: np.ndarray, seed: int, repeats: int = 10_000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(repeats, len(values)))
    means = values[indices].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def summarize_arm(arm: str, deltas: np.ndarray) -> dict[str, object]:
    confidence_low, confidence_high = bootstrap_mean_ci(deltas, seed=108215 + (0 if arm == "HEP" else 1))
    paired_t = stats.ttest_1samp(deltas, popmean=0.0)
    wilcoxon_p = np.nan if np.allclose(deltas, 0.0) else float(stats.wilcoxon(deltas, method="auto").pvalue)
    return {
        "analysis": "within_arm_paired_change",
        "arm": arm,
        "n_subjects": int(len(deltas)),
        "mean_post_minus_pre": float(deltas.mean()),
        "median_post_minus_pre": float(np.median(deltas)),
        "sd_post_minus_pre": float(np.std(deltas, ddof=1)),
        "mean_change_ci_low": confidence_low,
        "mean_change_ci_high": confidence_high,
        "paired_t_p_value_exploratory": float(paired_t.pvalue),
        "wilcoxon_p_value_exploratory": wilcoxon_p,
    }


def main() -> int:
    for source in [MATRIX, ANNOTATION, PANEL]:
        if not source.exists():
            raise FileNotFoundError(f"required source missing: {source}")

    panel = pd.read_csv(PANEL)
    panel["gene"] = panel["gene"].astype(str).str.strip().str.upper()
    if len(panel) != 12 or panel["gene"].nunique() != 12:
        raise RuntimeError("frozen panel must contain exactly 12 unique genes")
    if set(panel["trajectory_class"]) != {"MCI_peak", "MCI_trough"}:
        raise RuntimeError("frozen panel direction classes are malformed")

    metadata, probe_matrix = read_series_matrix()
    annotation = read_gpl570_annotation()
    genes = panel["gene"].tolist()
    probe_map = exact_probe_map(annotation, genes)
    mapped = [gene for gene in genes if probe_map[gene]]
    missing = [gene for gene in genes if not probe_map[gene]]
    if missing != EXPECTED_MISSING or len(mapped) != 10:
        raise RuntimeError(f"unexpected exact-symbol mapping: mapped={mapped}; missing={missing}")

    sample_columns = metadata["geo_accession"].tolist()
    expression = pd.DataFrame(index=sample_columns)
    for gene in mapped:
        values = probe_matrix.loc[probe_matrix["ID_REF"].astype(str).isin(probe_map[gene]), sample_columns]
        values = values.apply(pd.to_numeric, errors="coerce")
        if values.empty or values.isna().all(axis=None) or (values < 0).any(axis=None):
            raise RuntimeError(f"invalid GPL570 signal values for {gene}")
        expression[gene] = np.log2(values + 1.0).median(axis=0).to_numpy()
    if expression.isna().any().any():
        raise RuntimeError("partial score expression contains missing values")

    # All 86 rows are standardized before arm and timepoint labels are consulted.
    standardized = (expression - expression.mean(axis=0)) / expression.std(axis=0, ddof=1)
    directions = panel.set_index("gene")["trajectory_class"].to_dict()
    weights = pd.Series(
        {gene: (1.0 / 6.0 if directions[gene] == "MCI_peak" else -1.0 / 6.0) for gene in mapped}
    )
    score = standardized.loc[:, mapped].mul(weights, axis=1).sum(axis=1)

    metadata = metadata.copy()
    parsed = metadata["title"].map(parse_title)
    metadata[["arm", "timepoint", "subject"]] = pd.DataFrame(parsed.tolist(), index=metadata.index)
    metadata["score"] = metadata["geo_accession"].map(score)
    if metadata["score"].isna().any():
        raise RuntimeError("not every GEO sample received a frozen partial score")

    paired_deltas: dict[str, list[float]] = {"HEP": [], "MAP": []}
    pre_only: dict[str, int] = {"HEP": 0, "MAP": 0}
    post_only: dict[str, int] = {"HEP": 0, "MAP": 0}
    for (arm, _subject), rows in metadata.groupby(["arm", "subject"], sort=True):
        timepoints = set(rows["timepoint"])
        if len(rows) != len(timepoints):
            raise RuntimeError("GSE108215 contains duplicate observations for one arm/subject/timepoint")
        if timepoints == {"Pre", "Post"}:
            pre = float(rows.loc[rows["timepoint"] == "Pre", "score"].iloc[0])
            post = float(rows.loc[rows["timepoint"] == "Post", "score"].iloc[0])
            paired_deltas[arm].append(post - pre)
        elif timepoints == {"Pre"}:
            pre_only[arm] += 1
        elif timepoints == {"Post"}:
            post_only[arm] += 1
        else:
            raise RuntimeError("GSE108215 subject has an unsupported timepoint structure")

    if {arm: len(values) for arm, values in paired_deltas.items()} != {"HEP": 18, "MAP": 17}:
        raise RuntimeError("GSE108215 paired subject counts differ from the prespecified audit")
    hep = np.asarray(paired_deltas["HEP"], dtype=float)
    map_values = np.asarray(paired_deltas["MAP"], dtype=float)
    comparison = stats.ttest_ind(map_values, hep, equal_var=False)
    summaries = [summarize_arm("HEP", hep), summarize_arm("MAP", map_values)]
    summaries.append(
        {
            "analysis": "between_arm_difference_of_paired_changes",
            "arm": "MAP_minus_HEP",
            "n_subjects": int(len(map_values) + len(hep)),
            "mean_post_minus_pre": float(map_values.mean() - hep.mean()),
            "median_post_minus_pre": np.nan,
            "sd_post_minus_pre": np.nan,
            "mean_change_ci_low": np.nan,
            "mean_change_ci_high": np.nan,
            "paired_t_p_value_exploratory": float(comparison.pvalue),
            "wilcoxon_p_value_exploratory": np.nan,
        }
    )
    summary = pd.DataFrame(summaries)

    audit = {
        "status": "EXPLORATORY_PAIRED_INTERVENTION_RESPONSE_ONLY",
        "dataset": "GSE108215",
        "assay": "whole-blood Affymetrix U133 Plus 2.0 processed expression matrix",
        "study_design": "randomized MAP versus HEP MCI intervention cohort with pre/post blood samples",
        "analysis_unit": "paired subject-level post-minus-pre change",
        "sample_counts": {"total_expression_rows": 86, "pre": 50, "post": 36},
        "reported_vs_matrix_sample_counts": {
            "reported_pre": 51,
            "matrix_pre": 50,
            "reported_post": 36,
            "matrix_post": 36,
        },
        "sample_count_boundary": "The GEO series narrative reports 51 pre-treatment and 36 post-treatment samples; the official series matrix contains 50 pre-treatment and 36 post-treatment expression columns. The paired audit uses matrix-observed rows only.",
        "paired_subject_counts": {arm: len(values) for arm, values in paired_deltas.items()},
        "unpaired_subject_counts": {"pre_only": pre_only, "post_only": post_only},
        "mapping": {
            "panel_size": 12,
            "mapped": len(mapped),
            "missing_genes": missing,
            "probe_aggregation": "median of log2(signal+1) across exact-symbol GPL570 probes",
            "alias_substitution_allowed": False,
        },
        "score_definition": "frozen signed partial score: each mapped MCI_peak gene +1/6 and each mapped MCI_trough gene -1/6 after gene-wise standardization across all 86 samples",
        "arm_or_time_labels_used_for_preprocessing": False,
        "arm_or_time_labels_used_for_score_definition": False,
        "output_restriction": "aggregate results only; no GSM accessions, local participant tokens, individual scores, or expression values are exported",
        "source_hashes": {str(path): sha256(path) for path in [MATRIX, ANNOTATION, PANEL]},
        "official_sources": {
            "geo_record": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108215",
            "series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108215/matrix/GSE108215_series_matrix.txt.gz",
            "platform_annotation": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
            "access_date": "2026-08-13",
        },
        "allowed_claim": "This is an exploratory, label-leakage-free paired molecular response audit of a partial frozen score in a small public MCI intervention cohort.",
        "forbidden_claim": "This is not MCI-versus-HC validation, a complete 12-gene score replication, evidence of cognitive benefit, treatment efficacy, clinical utility, causality, or mechanism.",
    }

    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT / "gse108215_paired_summary.csv", index=False)
    summary.to_csv(SOURCE_OUT / "gse108215_paired_summary.csv", index=False)
    for destination in [OUT / "gse108215_paired_audit.json", SOURCE_OUT / "gse108215_paired_audit.json"]:
        destination.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
