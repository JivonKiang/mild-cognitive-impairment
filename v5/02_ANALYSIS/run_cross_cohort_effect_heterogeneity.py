from __future__ import annotations

import argparse
import json
import re
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from run_cross_study_stable_signal import development_endpoint, load_gse140829_unlabelled
from run_gse249477_rank_transport import (
    GSE249477_RAW,
    GSE249477_SERIES,
    parse_assessment_labels as parse_gse249477_assessment_labels,
)
from run_pairwise_external import (
    GPL15988_FULL,
    GSE140829_MATRIX,
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    load_external,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
LABELS = {"HC": 0, "MCI": 1}
EXPECTED_SHARED_GENE_COUNT = 11575
COHORT_ORDER = ["GSE63060", "GSE63061", "GSE140829", "GSE249477"]


def load_gse249477_unlabelled_expression() -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    """Map GSE249477 TPM columns without opening diagnosis-bearing series titles."""
    raw = pd.read_csv(GSE249477_RAW, compression="gzip", low_memory=False)
    tpm_columns = [column for column in raw.columns if " - TPM" in column]
    if len(tpm_columns) != 62:
        raise ValueError(f"expected 62 GSE249477 TPM columns, found {len(tpm_columns)}")
    sample_ids: list[str] = []
    for column in tpm_columns:
        match = re.fullmatch(r"(DK22119_\d+) \(GE\) - TPM", column)
        if match is None:
            raise ValueError(f"TPM column has unexpected GSE249477 header: {column!r}")
        sample_ids.append(match.group(1))
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("GSE249477 TPM identifiers are not unique")
    values = raw.loc[:, tpm_columns].apply(pd.to_numeric, errors="coerce")
    values.columns = sample_ids
    symbols = raw["Name"].astype(str).str.upper()
    expression = values.groupby(symbols, sort=True).median().T.loc[sample_ids]
    duplicate_rows_collapsed = int(len(raw) - expression.shape[1])
    if duplicate_rows_collapsed != 1:
        raise ValueError(f"expected exactly one duplicate target gene-row collapse, found {duplicate_rows_collapsed}")
    return expression.astype(np.float32), sample_ids, {
        "target_tpm_columns": len(sample_ids),
        "unique_tpm_ids": len(set(sample_ids)),
        "target_gene_count_after_exact_symbol_median_collapse": int(expression.shape[1]),
        "target_duplicate_gene_rows_collapsed": duplicate_rows_collapsed,
        "series_titles_read_before_shared_universe": False,
    }


def validate_effect_universe(effect_series: dict[str, pd.Series]) -> list[str]:
    """Require every cohort effect vector to have the same named gene set."""
    if not effect_series:
        raise ValueError("at least one cohort effect vector is required")
    reference: set[str] | None = None
    for cohort, values in effect_series.items():
        if not values.index.is_unique or not all(isinstance(gene, str) for gene in values.index):
            raise ValueError("all cohort effect vectors must use the same gene universe with named unique genes")
        current = set(values.index)
        if not current:
            raise ValueError("all cohort effect vectors must use the same gene universe with named unique genes")
        if reference is None:
            reference = current
        elif current != reference:
            raise ValueError("all cohort effect vectors must use the same gene universe")
    return sorted(reference or [])


def cohort_effects(expression: pd.DataFrame, y: np.ndarray, cohort: str) -> pd.DataFrame:
    """Compute independent Hedges' g (MCI minus HC) and its sampling variance."""
    if len(expression) != len(y):
        raise ValueError("expression rows and endpoint labels must have the same length")
    if not expression.columns.is_unique or not all(isinstance(gene, str) for gene in expression.columns):
        raise ValueError("expression requires unique named gene columns")
    values = expression.to_numpy(dtype=float)
    cases = values[np.asarray(y, dtype=int) == 1]
    controls = values[np.asarray(y, dtype=int) == 0]
    n_mci, n_hc = len(cases), len(controls)
    if n_mci < 2 or n_hc < 2:
        raise ValueError("each cohort requires at least two MCI and two HC rows")
    def group_moments(group: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        finite = np.isfinite(group)
        count = finite.sum(axis=0)
        total = np.where(finite, group, 0.0).sum(axis=0)
        mean = np.divide(total, count, out=np.full(group.shape[1], np.nan), where=count > 0)
        squared_deviation = np.where(finite, (group - mean) ** 2, 0.0).sum(axis=0)
        variance = np.divide(
            squared_deviation,
            count - 1,
            out=np.full(group.shape[1], np.nan),
            where=count > 1,
        )
        return count, mean, variance

    mci_count, mci_mean, mci_variance = group_moments(cases)
    hc_count, hc_mean, hc_variance = group_moments(controls)
    degrees_of_freedom = mci_count + hc_count - 2
    correction = np.divide(
        3.0,
        4.0 * degrees_of_freedom - 1.0,
        out=np.full(expression.shape[1], np.nan),
        where=degrees_of_freedom > 0,
    )
    correction = 1.0 - correction
    pooled_standard_deviation = np.sqrt(
        np.divide(
            (mci_count - 1) * mci_variance + (hc_count - 1) * hc_variance,
            degrees_of_freedom,
            out=np.full(expression.shape[1], np.nan),
            where=degrees_of_freedom > 0,
        )
    )
    standardized_difference = np.divide(
        mci_mean - hc_mean,
        pooled_standard_deviation,
        out=np.full(expression.shape[1], np.nan, dtype=float),
        where=np.isfinite(pooled_standard_deviation) & (pooled_standard_deviation > 0),
    )
    hedges_g = correction * standardized_difference
    total_count = mci_count + hc_count
    sampling_variance = np.divide(
        total_count,
        mci_count * hc_count,
        out=np.full(expression.shape[1], np.nan),
        where=(mci_count > 0) & (hc_count > 0),
    ) + np.divide(
        hedges_g**2,
        2.0 * total_count,
        out=np.full(expression.shape[1], np.nan),
        where=total_count > 0,
    )
    estimable = (
        (mci_count >= 2)
        & (hc_count >= 2)
        & np.isfinite(hedges_g)
        & np.isfinite(sampling_variance)
        & (sampling_variance > 0)
    )
    result = pd.DataFrame(
        {
            "cohort": cohort,
            "gene": expression.columns,
            "n_mci": mci_count,
            "n_hc": hc_count,
            "hedges_g_mci_minus_hc": hedges_g,
            "sampling_variance": sampling_variance,
            "effect_sign": np.where(hedges_g > 0, "MCI_higher", np.where(hedges_g < 0, "HC_higher", "zero_or_not_estimable")),
            "status": np.where(estimable, "estimable", "non_estimable_zero_or_nonfinite_variance"),
        }
    )
    return result


def pairwise_effect_summary(effects: pd.DataFrame, cohort_a: str, cohort_b: str) -> dict[str, Any]:
    """Summarize effect concordance without fitting a cross-cohort model."""
    required = {"gene", cohort_a, cohort_b}
    if not required.issubset(effects.columns):
        raise ValueError(f"effects requires columns {sorted(required)}")
    if effects["gene"].duplicated().any() or not all(isinstance(gene, str) for gene in effects["gene"]):
        raise ValueError("effects requires unique named genes")
    values = effects[[cohort_a, cohort_b]].apply(pd.to_numeric, errors="coerce")
    usable = np.isfinite(values[cohort_a]) & np.isfinite(values[cohort_b])
    left = values.loc[usable, cohort_a].to_numpy(dtype=float)
    right = values.loc[usable, cohort_b].to_numpy(dtype=float)
    if len(left) < 2 or np.ptp(left) == 0 or np.ptp(right) == 0:
        correlation = np.nan
    else:
        correlation = float(spearmanr(left, right).statistic)
    nonzero = (left != 0.0) & (right != 0.0)
    sign_concordance = float(np.mean(np.sign(left[nonzero]) == np.sign(right[nonzero]))) if nonzero.any() else np.nan
    return {
        "cohort_a": cohort_a,
        "cohort_b": cohort_b,
        "comparison_type": "cohort_pair_descriptive",
        "estimable_gene_count": int(len(left)),
        "nonzero_both_gene_count": int(nonzero.sum()),
        "spearman_effect_correlation": correlation,
        "sign_concordance_fraction": sign_concordance,
        "inference_boundary": "descriptive_only_no_pooled_effect",
    }


def development_heterogeneity(effects: pd.DataFrame) -> pd.DataFrame:
    """Two-study AddNeuroMed Q and I-squared, with no released pooled estimate."""
    required = {
        "gene",
        "GSE63060_effect",
        "GSE63061_effect",
        "GSE63060_variance",
        "GSE63061_variance",
    }
    if not required.issubset(effects.columns):
        raise ValueError("development heterogeneity input is incomplete")
    left_effect = effects["GSE63060_effect"].to_numpy(dtype=float)
    right_effect = effects["GSE63061_effect"].to_numpy(dtype=float)
    left_variance = effects["GSE63060_variance"].to_numpy(dtype=float)
    right_variance = effects["GSE63061_variance"].to_numpy(dtype=float)
    valid = (
        np.isfinite(left_effect)
        & np.isfinite(right_effect)
        & np.isfinite(left_variance)
        & np.isfinite(right_variance)
        & (left_variance > 0)
        & (right_variance > 0)
    )
    left_weight = np.divide(1.0, left_variance, out=np.full(len(effects), np.nan), where=valid)
    right_weight = np.divide(1.0, right_variance, out=np.full(len(effects), np.nan), where=valid)
    internal_reference = np.divide(
        left_weight * left_effect + right_weight * right_effect,
        left_weight + right_weight,
        out=np.full(len(effects), np.nan),
        where=valid,
    )
    q = left_weight * (left_effect - internal_reference) ** 2 + right_weight * (right_effect - internal_reference) ** 2
    i_squared = np.full(len(effects), np.nan)
    positive_q = valid & (q > 0)
    i_squared[positive_q] = np.maximum(0.0, (q[positive_q] - 1.0) / q[positive_q]) * 100.0
    return pd.DataFrame(
        {
            "gene": effects["gene"],
            "study_count": 2,
            "cochran_q": np.where(valid, q, np.nan),
            "i_squared_percent": i_squared,
            "status": np.where(valid, "estimable_two_study_descriptive", "non_estimable"),
            "inference_boundary": "AddNeuroMed_only_no_released_pooled_effect_or_p_value",
        }
    )


def _endpoint_expression() -> tuple[dict[str, pd.DataFrame], dict[str, np.ndarray], dict[str, Any]]:
    x60, y60, x61, y61, development_genes = development_endpoint()
    development_universe = set(development_genes)

    # Both external expression matrices are mapped before any target diagnosis is read.
    x140_unlabelled, audit140_unlabelled = load_gse140829_unlabelled(development_universe)
    x249_all, tpm_ids, audit249_unlabelled = load_gse249477_unlabelled_expression()
    shared_genes = sorted(development_universe.intersection(x140_unlabelled.columns, x249_all.columns))
    if len(shared_genes) != EXPECTED_SHARED_GENE_COUNT:
        raise ValueError(f"expected {EXPECTED_SHARED_GENE_COUNT} exact shared genes, found {len(shared_genes)}")

    x140_labelled, meta140, audit140_labelled = load_external(
        GSE140829_MATRIX, GPL15988_FULL, set(shared_genes)
    )
    keep140 = meta140["label"].isin(LABELS)
    x140 = x140_labelled.loc[keep140, shared_genes]
    y140 = meta140.loc[x140.index, "label"].map(LABELS).to_numpy(dtype=int)
    if (len(x140), int(y140.sum())) != (383, 134):
        raise ValueError("unexpected GSE140829 MCI/HC sample-row endpoint counts")

    labels249 = parse_gse249477_assessment_labels(tpm_ids)
    keep249 = labels249.isin(LABELS)
    x249 = x249_all.loc[keep249, shared_genes]
    y249 = labels249.loc[keep249].map(LABELS).to_numpy(dtype=int)
    if (len(x249), int(y249.sum())) != (41, 20):
        raise ValueError("unexpected GSE249477 MCI/HC endpoint counts")

    expression = {
        "GSE63060": x60.loc[:, shared_genes],
        "GSE63061": x61.loc[:, shared_genes],
        "GSE140829": x140,
        "GSE249477": x249,
    }
    endpoint = {"GSE63060": y60, "GSE63061": y61, "GSE140829": y140, "GSE249477": y249}
    audit = {
        "shared_exact_case_canonicalized_gene_count": len(shared_genes),
        "shared_universe_built_before_external_diagnosis_labels": True,
        "gse140829_unlabelled_expression_audit": audit140_unlabelled,
        "gse140829_assessment_metadata_audit": audit140_labelled,
        "gse249477_unlabelled_expression_audit": audit249_unlabelled,
    }
    return expression, endpoint, audit


def main(output_dir: Path = OUT) -> int:
    expression, endpoint, input_audit = _endpoint_expression()
    per_cohort = [cohort_effects(expression[cohort], endpoint[cohort], cohort) for cohort in COHORT_ORDER]
    effect_by_cohort = {
        cohort: table.set_index("gene")["hedges_g_mci_minus_hc"]
        for cohort, table in zip(COHORT_ORDER, per_cohort)
    }
    shared_genes = validate_effect_universe(effect_by_cohort)
    all_effects = pd.concat(per_cohort, ignore_index=True).sort_values(["gene", "cohort"], kind="mergesort")
    all_effects["cohort"] = pd.Categorical(all_effects["cohort"], categories=COHORT_ORDER, ordered=True)
    all_effects = all_effects.sort_values(["gene", "cohort"], kind="mergesort")
    all_effects["cohort"] = all_effects["cohort"].astype(str)

    wide_effects = pd.DataFrame({"gene": shared_genes})
    for cohort, table in zip(COHORT_ORDER, per_cohort):
        indexed = table.set_index("gene")
        wide_effects[cohort] = indexed.loc[shared_genes, "hedges_g_mci_minus_hc"].to_numpy()
        if cohort in {"GSE63060", "GSE63061"}:
            wide_effects[f"{cohort}_effect"] = indexed.loc[shared_genes, "hedges_g_mci_minus_hc"].to_numpy()
            wide_effects[f"{cohort}_variance"] = indexed.loc[shared_genes, "sampling_variance"].to_numpy()

    pairwise = [pairwise_effect_summary(wide_effects, left, right) for left, right in combinations(COHORT_ORDER, 2)]
    development_reference = (wide_effects["GSE63060"] + wide_effects["GSE63061"]) / 2.0
    reference_effects = pd.DataFrame({"gene": shared_genes, "AddNeuroMed_unweighted_mean": development_reference})
    for cohort in ["GSE140829", "GSE249477"]:
        summary = pairwise_effect_summary(
            reference_effects.assign(**{cohort: wide_effects[cohort]}), "AddNeuroMed_unweighted_mean", cohort
        )
        summary["comparison_type"] = "development_family_unweighted_reference_descriptive"
        summary["inference_boundary"] = "descriptive_reference_no_pooling_or_cross_cohort_inference"
        pairwise.append(summary)
    heterogeneity = development_heterogeneity(wide_effects)

    output_dir.mkdir(parents=True, exist_ok=True)
    all_effects.to_csv(output_dir / "cross_cohort_effects.csv", index=False)
    pd.DataFrame(pairwise).to_csv(output_dir / "cross_cohort_effect_pairwise_summary.csv", index=False)
    heterogeneity.to_csv(output_dir / "cross_cohort_effect_heterogeneity.csv", index=False)
    audit = {
        "status": "CROSS_COHORT_EFFECT_HETEROGENEITY_AUDIT",
        "endpoint": "MCI_minus_HC",
        "analysis_type": "nonpredictive_descriptive_effect_concordance",
        "cohort_order": COHORT_ORDER,
        "cohort_endpoint_counts": {
            cohort: {"sample_rows": int(len(endpoint[cohort])), "mci_rows": int(endpoint[cohort].sum()), "hc_rows": int((endpoint[cohort] == 0).sum())}
            for cohort in COHORT_ORDER
        },
        **input_audit,
        "gse140829_subject_independence_resolved": False,
        "gse140829_excluded_from_pooled_inference": True,
        "external_labels_used_for_mapping": False,
        "external_labels_used_for_feature_selection": False,
        "external_labels_used_for_normalization": False,
        "external_labels_used_for_model_fitting": False,
        "external_labels_used_for_calibration_thresholding_or_adaptation": False,
        "external_labels_used_only_for_independent_cohort_effect_assessment": True,
        "development_heterogeneity_scope": "GSE63060_and_GSE63061_only; two-study descriptive Q and I-squared",
        "no_pooled_four_cohort_effect": True,
        "no_cross_cohort_heterogeneity_p_value": True,
        "allowed_claim": "This audit describes observed gene-effect concordance or discordance in the locked shared gene universe.",
        "forbidden_claim": "This audit proves a biological mechanism, blood-cell causality, a clinically useful biomarker, or successful transportability.",
        "method_reference": {
            "citation": "Higgins et al., BMJ 2003;327:557",
            "pmid": "12958120",
            "doi": "10.1136/bmj.327.7414.557",
            "url": "https://pubmed.ncbi.nlm.nih.gov/12958120/",
            "accessed": "2026-08-13",
        },
        "official_source_urls": {
            "GSE140829": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829",
            "GSE249477": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249477",
        },
        "source_hashes": {
            "GSE63060_data_matrix.csv": sha256_file(GSE60_MATRIX),
            "GSE63060_official_sample_labels.csv": sha256_file(GSE60_LABELS),
            "GSE63061_data_matrix.csv": sha256_file(GSE61_MATRIX),
            "GSE63061_official_sample_labels.csv": sha256_file(GSE61_LABELS),
            "GSE140829_series_matrix.fresh.txt.gz": sha256_file(GSE140829_MATRIX),
            "GPL15988_full.txt": sha256_file(GPL15988_FULL),
            "GSE249477_raw_count_normalize_04-10-2025.csv.gz": sha256_file(GSE249477_RAW),
            "GSE249477_series_matrix.txt.gz": sha256_file(GSE249477_SERIES),
            "run_cross_cohort_effect_heterogeneity.py": sha256_file(Path(__file__).resolve()),
        },
    }
    (output_dir / "cross_cohort_effect_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"shared_genes": len(shared_genes), "effects": len(all_effects), "pairwise_summaries": len(pairwise)}, indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUT)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
