from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from run_pairwise_external import (
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    bootstrap_auc,
    load_dev,
    metric_row,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260813
GSE249477_RAW = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_raw_count_normalize_04-10-2025.csv.gz"
)
GSE249477_SERIES = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_series_matrix.txt.gz"
)


def canonical_symbol(value: object) -> str:
    return str(value).upper()


def sample_rank_transform(frame: pd.DataFrame) -> pd.DataFrame:
    """Rank and z-score each sample independently, without target labels."""
    ranks = np.apply_along_axis(
        lambda row: pd.Series(row).rank(method="average").to_numpy(), 1, frame.to_numpy(dtype=float)
    )
    center = np.nanmean(ranks, axis=1, keepdims=True)
    scale = np.nanstd(ranks, axis=1, keepdims=True)
    return pd.DataFrame(
        (ranks - center) / np.where(scale > 0, scale, 1.0),
        index=frame.index,
        columns=frame.columns,
    ).astype(np.float32)


def parse_series_title_ids() -> list[str]:
    titles: list[str] = []
    with gzip.open(GSE249477_SERIES, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:]
                break
    if len(titles) != 62:
        raise ValueError(f"expected 62 GSE249477 series titles, found {len(titles)}")
    identifiers: list[str] = []
    for title in titles:
        matches = re.findall(r"\[\s*(DK22119_\d+)\s*\]", title)
        if len(matches) != 1:
            raise ValueError(f"title must contain exactly one GSE249477 identifier: {title!r}")
        identifiers.append(matches[0])
    if len(set(identifiers)) != 62:
        raise ValueError("GSE249477 series title identifiers are not unique")
    return identifiers


def parse_assessment_labels(series_ids: list[str]) -> pd.Series:
    """Access GSE249477 diagnosis labels only after predictions are already fixed."""
    titles: list[str] = []
    with gzip.open(GSE249477_SERIES, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:]
                break
    labels: dict[str, str] = {}
    for title in titles:
        matches = re.findall(r"\[\s*(DK22119_\d+)\s*\]", title)
        if len(matches) != 1:
            raise ValueError(f"title must contain exactly one GSE249477 identifier: {title!r}")
        sample_id = matches[0]
        lower = title.casefold()
        if "mild cognitive impairment" in lower:
            labels[sample_id] = "MCI"
        elif "cognitively normal" in lower:
            labels[sample_id] = "HC"
        elif "alzheimer" in lower:
            labels[sample_id] = "AD"
        else:
            raise ValueError(f"title lacks recognized diagnosis: {title!r}")
    parsed = pd.Series(labels, name="label").reindex(series_ids)
    if parsed.isna().any():
        raise ValueError("missing GSE249477 diagnosis after title parsing")
    return parsed


def load_target_expression() -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    series_ids = parse_series_title_ids()
    raw = pd.read_csv(GSE249477_RAW, compression="gzip", low_memory=False)
    tpm_columns = [column for column in raw.columns if " - TPM" in column]
    if len(tpm_columns) != 62:
        raise ValueError(f"expected 62 TPM columns, found {len(tpm_columns)}")
    column_ids: list[str] = []
    for column in tpm_columns:
        match = re.fullmatch(r"(DK22119_\d+) \(GE\) - TPM", column)
        if match is None:
            raise ValueError(f"TPM column has unexpected GSE249477 header: {column!r}")
        column_ids.append(match.group(1))
    if len(set(column_ids)) != 62:
        raise ValueError("GSE249477 TPM identifiers are not unique")
    if set(column_ids) != set(series_ids):
        raise ValueError("GSE249477 title and TPM identifiers do not match")

    values = raw.loc[:, tpm_columns].apply(pd.to_numeric, errors="coerce")
    values.columns = column_ids
    raw["symbol"] = raw["Name"].map(canonical_symbol)
    expression = values.groupby(raw["symbol"], sort=True).median().T
    expression = expression.loc[column_ids]
    duplicate_rows_collapsed = int(len(raw) - expression.shape[1])
    if duplicate_rows_collapsed != 1:
        raise ValueError(f"expected exactly one duplicate target gene-row collapse, found {duplicate_rows_collapsed}")
    audit = {
        "series_title_count": int(len(series_ids)),
        "target_tpm_columns": int(len(column_ids)),
        "unique_title_ids": int(len(set(series_ids))),
        "unique_tpm_ids": int(len(set(column_ids))),
        "exact_id_mapping": True,
        "target_gene_count_after_exact_symbol_median_collapse": int(expression.shape[1]),
        "target_duplicate_gene_rows_collapsed": duplicate_rows_collapsed,
    }
    return expression.astype(np.float32), series_ids, audit


def make_estimator(feature_count: int, seed: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("select", SelectKBest(score_func=f_classif, k=min(200, feature_count))),
            ("scale", StandardScaler()),
            ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
        ]
    )


def load_development() -> tuple[pd.DataFrame, pd.DataFrame]:
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    common = sorted(set(dev60_x.columns).intersection(dev61_x.columns))
    expression = pd.concat([dev60_x[common], dev61_x[common]], axis=0)
    metadata = pd.concat([dev60_meta, dev61_meta], axis=0).loc[expression.index]
    keep = metadata["label"].isin({"HC", "MCI"})
    return expression.loc[keep], metadata.loc[keep]


def generate_target_predictions(target_expression: pd.DataFrame, model: Pipeline) -> tuple[pd.DataFrame, np.ndarray]:
    """Return all-target rank/z features and predictions without labels or cohort fitting."""
    transformed = sample_rank_transform(target_expression)
    return transformed, model.predict_proba(transformed)[:, 1]


def target_assessment(y: np.ndarray, probability: np.ndarray) -> dict[str, float | int | str]:
    """Compute final target metrics without fitting target-label calibration or thresholds."""
    if len(np.unique(y)) != 2:
        raise ValueError("target assessment requires both HC and MCI")
    bins = np.linspace(0.0, 1.0, 11)
    ece = 0.0
    used_bins = 0
    for low, high in zip(bins[:-1], bins[1:]):
        mask = (probability >= low) & (probability < high if high < 1 else probability <= high)
        if mask.any():
            ece += float(mask.mean()) * abs(float(probability[mask].mean()) - float(y[mask].mean()))
            used_bins += 1
    ci_low, ci_high = bootstrap_auc(y, probability, SEED + 200)
    return {
        "status": "SECONDARY_INDEPENDENT_TRANSPORT_AUDIT",
        "endpoint": "stage_mci_vs_hc",
        "cohort": "GSE249477",
        "model": "sample_rank_transform",
        "n": int(len(y)),
        "positive_n": int(y.sum()),
        "auc": float(roc_auc_score(y, probability)),
        "auc_ci_low": ci_low,
        "auc_ci_high": ci_high,
        "brier": float(brier_score_loss(y, probability)),
        "ece": float(ece),
        "ece_nonempty_bins": used_bins,
        "calibration_method": "NOT_FIT_TARGET_LABELS",
        "calibration_fit_scope": "not_fit_target_labels",
    }


def main(output_dir: Path = OUT) -> int:
    dev_expression, dev_metadata = load_development()
    target_expression, _series_ids, target_input_audit = load_target_expression()
    common_genes = sorted(set(dev_expression.columns).intersection(target_expression.columns))
    if not common_genes:
        raise ValueError("no exact case-canonicalized genes overlap development and GSE249477")

    dev_ranked = sample_rank_transform(dev_expression[common_genes])
    y_dev = dev_metadata["label"].map({"HC": 0, "MCI": 1}).to_numpy(dtype=int)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = np.full(len(y_dev), np.nan, dtype=float)
    for fold, (train_idx, test_idx) in enumerate(cv.split(dev_ranked, y_dev), start=1):
        model = make_estimator(dev_ranked.shape[1], SEED + fold)
        model.fit(dev_ranked.iloc[train_idx], y_dev[train_idx])
        oof[test_idx] = model.predict_proba(dev_ranked.iloc[test_idx])[:, 1]
    final_model = make_estimator(dev_ranked.shape[1], SEED + 999)
    final_model.fit(dev_ranked, y_dev)

    # All 62 target rows are transformed and predicted before any diagnosis label is parsed.
    target_ranked_all, target_probability_all = generate_target_predictions(target_expression[common_genes], final_model)

    # Target labels are accessed only to lock and quantify the final HC/MCI assessment.
    target_labels = parse_assessment_labels(target_ranked_all.index.tolist())
    assessment_mask = target_labels.isin({"HC", "MCI"})
    assessment_labels = target_labels.loc[assessment_mask]
    y_target = assessment_labels.map({"HC": 0, "MCI": 1}).to_numpy(dtype=int)
    target_probability = target_probability_all[assessment_mask.to_numpy()]

    internal = metric_row(
        y_dev,
        oof,
        "NESTED_INTERNAL_SENSITIVITY",
        "stage_mci_vs_hc",
        "GSE63060_GSE63061",
        "sample_rank_transform",
    )
    internal["auc_ci_low"], internal["auc_ci_high"] = bootstrap_auc(y_dev, oof, SEED + 100)
    internal["transformation"] = "within-sample rank then within-sample z-score"
    external = target_assessment(y_target, target_probability)
    external["transformation"] = "within-sample rank then within-sample z-score"
    external["target_labels_used_for_transform"] = False
    external["target_labels_used_for_fit"] = False
    external["target_labels_used_for_feature_selection"] = False
    external["target_labels_used_for_threshold_or_calibration"] = False
    external["target_labels_used_only_for_assessment"] = True

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([internal, external]).to_csv(output_dir / "gse249477_rank_transport_metrics.csv", index=False)
    pd.DataFrame(
        {
            "cohort": "GSE249477",
            "assessment_index": np.arange(len(y_target), dtype=int),
            "y": y_target,
            "probability": target_probability,
        }
    ).to_csv(output_dir / "gse249477_rank_transport_predictions_anonymized.csv", index=False)
    audit = {
        "status": "SECONDARY_INDEPENDENT_TRANSPORT_AUDIT",
        "assessment_cohort": "GSE249477",
        "assessment_unit": "sample-row; 21 HC and 20 MCI",
        "assessment_rows": int(len(assessment_labels)),
        "assessment_label_counts": dict(sorted(assessment_labels.value_counts().to_dict().items())),
        "excluded_ad_rows": int((target_labels == "AD").sum()),
        "transformed_target_rows": int(len(target_ranked_all)),
        "target_labels_used_for_transform": False,
        "target_labels_used_for_fit": False,
        "target_labels_used_for_feature_selection": False,
        "target_labels_used_for_threshold_or_calibration": False,
        "target_labels_used_only_for_assessment": True,
        **target_input_audit,
        "target_identifier_audit": target_input_audit,
        "development": {
            "datasets": ["GSE63060", "GSE63061"],
            "row_count": int(len(dev_metadata)),
            "label_counts": dict(sorted(dev_metadata["label"].value_counts().to_dict().items())),
        },
        "common_exact_case_canonicalized_gene_count": int(len(common_genes)),
        "model_spec": {
            "development_cv": "5-fold StratifiedKFold; AddNeuroMed only",
            "final_model_fit_scope": "GSE63060 and GSE63061 HC/MCI rows only",
            "transformation": "within-sample rank then within-sample z-score",
            "pipeline": "SimpleImputer(median), SelectKBest(f_classif, k=min(200,n_features)), StandardScaler, LogisticRegression(C=0.2, ridge)",
        },
        "allowed_claims": [
            "Secondary independent cross-platform transport audit of a prespecified rank/z-score sensitivity.",
            "GSE249477 diagnosis labels were used only to quantify the locked HC/MCI assessment.",
        ],
        "forbidden_claims": [
            "GSE249477 clinically validates the model.",
            "The rank transform solves domain shift.",
            "The results prove a generalizable MCI blood biomarker.",
        ],
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [
                GSE60_MATRIX,
                GSE60_LABELS,
                GSE61_MATRIX,
                GSE61_LABELS,
                GSE249477_RAW,
                GSE249477_SERIES,
                Path(__file__).resolve(),
            ]
        },
        "official_source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249477",
        "official_source_url_verification": "PENDING_ONLINE_VERIFICATION",
    }
    (output_dir / "gse249477_rank_transport_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(pd.DataFrame([internal, external]).to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUT)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
