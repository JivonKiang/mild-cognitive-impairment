from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
PROTOCOL = ROOT / "00_PROTOCOL" / "analysis_protocol_pairwise_v1.md"
PANEL = [
    ("GPI", 1),
    ("PHF15", 1),
    ("VAT1", 1),
    ("EIF3B", 1),
    ("UCP2", 1),
    ("LAT1-3TM", 1),
    ("C22ORF39", -1),
    ("USP38", -1),
    ("ZSWIM6", -1),
    ("PRDM10", -1),
    ("RP2", -1),
    ("GCA", -1),
]

GSE60_MATRIX = Path(r"E:\20241004_MCI\20250522\GEO\data\GSE63060_data_matrix.csv")
GSE60_LABELS = Path(r"E:\20241004_MCI\20260709\data\GSE63060_official_sample_labels.csv")
GSE61_MATRIX = Path(r"E:\20241004_MCI\20260709\data\GSE63061_data_matrix.csv")
GSE61_LABELS = Path(r"E:\20241004_MCI\20260709\data\GSE63061_official_sample_labels.csv")
GSE140829_MATRIX = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse140829_task8\GSE140829_series_matrix.fresh.txt.gz"
)
GPL15988_FULL = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse140829_task9\GPL15988_full.txt"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clean(value: str) -> str:
    return value.strip().strip('"')


def _canonical_symbol(value: str) -> str:
    return _clean(str(value)).upper()


def _label(value: str) -> str:
    value = str(value).strip().casefold()
    if value in {"ctl", "control", "hc", "normal"}:
        return "HC"
    if value == "mci" or "mild cognitive" in value:
        return "MCI"
    if value == "ad" or "alzheimer" in value:
        return "AD"
    return "OTHER"


def _numeric(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def load_dev(
    matrix_path: Path,
    labels_path: Path,
    dataset_id: str,
    include_unlabelled_rows: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(matrix_path)
    feature_col = "Name" if "Name" in matrix.columns else matrix.columns[0]
    feature_names = matrix[feature_col].astype(str).map(_canonical_symbol)
    matrix = matrix.drop(columns=[feature_col])
    matrix.index = feature_names
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.groupby(level=0).mean()
    expression = matrix.T

    labels = pd.read_csv(labels_path)
    labels["sample"] = labels["sample"].astype(str)
    included = labels["included"].astype(str).str.casefold().isin({"yes", "true", "1"})
    labels = labels.loc[included].copy()
    labels["label"] = labels["status"].map(_label)
    labels["age"] = pd.to_numeric(labels.get("age"), errors="coerce")
    labels["sex_num"] = labels.get("gender", "").astype(str).str.casefold().map(
        {"male": 0.0, "m": 0.0, "female": 1.0, "f": 1.0}
    )
    if not include_unlabelled_rows:
        labels = labels[labels["label"].isin({"HC", "MCI", "AD"})].copy()
    labels = labels[labels["sample"].isin(expression.index)].copy()
    labels = labels.drop_duplicates("sample", keep="first").set_index("sample")
    expression = expression.loc[labels.index]
    labels["dataset"] = dataset_id
    labels["study_family"] = "AddNeuroMed"
    expression = expression.astype(np.float32)
    return expression, labels[["label", "age", "sex_num", "dataset", "study_family"]]


def load_gpl_symbol_map(path: Path) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        in_table = False
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!platform_table_begin":
                next(handle)
                in_table = True
                continue
            if in_table and line == "!platform_table_end":
                break
            if not in_table:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            probe = _clean(fields[0])
            symbols = {
                _canonical_symbol(item)
                for item in fields[2].split("///")
                if re.fullmatch(r"[A-Za-z][A-Za-z0-9-]*", item.strip())
            }
            if probe and symbols:
                mapping[probe] = symbols
    return mapping


def _parse_characteristics(rows: list[list[str]], sample_count: int) -> tuple[list[str], np.ndarray, np.ndarray]:
    labels = ["OTHER"] * sample_count
    ages = np.full(sample_count, np.nan, dtype=float)
    sexes = np.full(sample_count, np.nan, dtype=float)
    for row in rows:
        for idx, raw in enumerate(row[:sample_count]):
            value = _clean(raw)
            lower = value.casefold()
            if "diagnosis:" in lower:
                labels[idx] = _label(value.split(":", 1)[1].strip())
            elif re.search(r"\b(mci|ad|control)\b", lower):
                labels[idx] = _label(value)
            age_match = re.search(r"age[^:]*:\s*([0-9]+(?:\.[0-9]+)?)", lower)
            if age_match:
                ages[idx] = float(age_match.group(1))
            if re.search(r"sex\s*:\s*(male|m)", lower):
                sexes[idx] = 0.0
            elif re.search(r"sex\s*:\s*(female|f)", lower):
                sexes[idx] = 1.0
    return labels, ages, sexes


def load_external_expression_unlabeled(
    series_path: Path,
    platform_path: Path,
    allowed_genes: set[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load all expression rows without reading sample diagnosis metadata."""
    probe_symbols = load_gpl_symbol_map(platform_path)
    sample_ids: list[str] = []
    gene_sum: dict[str, np.ndarray] = {}
    gene_count: dict[str, np.ndarray] = {}
    in_table = False
    with gzip.open(series_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not in_table:
                if line == "!series_matrix_table_begin":
                    header = next(csv.reader([next(handle)], delimiter="\t"))
                    sample_ids = [_clean(value) for value in header[1:]]
                    in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            fields = line.split("\t")
            probe = _clean(fields[0]) if fields else ""
            symbols = probe_symbols.get(probe, set()).intersection(allowed_genes)
            if not symbols:
                continue
            values = np.array([_numeric(value) for value in fields[1 : len(sample_ids) + 1]], dtype=np.float32)
            if values.size != len(sample_ids):
                continue
            valid = np.isfinite(values)
            for gene in symbols:
                if gene not in gene_sum:
                    gene_sum[gene] = np.zeros(len(sample_ids), dtype=np.float64)
                    gene_count[gene] = np.zeros(len(sample_ids), dtype=np.float64)
                gene_sum[gene][valid] += values[valid]
                gene_count[gene][valid] += 1.0

    expression = pd.DataFrame(
        {
            gene: np.divide(
                gene_sum[gene],
                gene_count[gene],
                out=np.full(len(sample_ids), np.nan, dtype=np.float64),
                where=gene_count[gene] > 0,
            )
            for gene in sorted(gene_sum)
        },
        index=sample_ids,
    ).astype(np.float32)
    audit = {
        "series_sample_count": len(sample_ids),
        "platform_probe_count": len(probe_symbols),
        "mapped_gene_count": len(gene_sum),
        "diagnosis_metadata_accessed": False,
    }
    return expression, audit


def load_external(
    series_path: Path,
    platform_path: Path,
    allowed_genes: set[str],
    include_unlabelled_rows: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    probe_symbols = load_gpl_symbol_map(platform_path)
    metadata: dict[str, list[list[str]]] = {}
    sample_ids: list[str] = []
    gene_sum: dict[str, np.ndarray] = {}
    gene_count: dict[str, np.ndarray] = {}
    in_table = False
    with gzip.open(series_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not in_table:
                if line == "!series_matrix_table_begin":
                    header = next(csv.reader([next(handle)], delimiter="\t"))
                    sample_ids = [_clean(value) for value in header[1:]]
                    in_table = True
                    continue
                fields = next(csv.reader([line], delimiter="\t"))
                if fields and fields[0].startswith("!Sample_"):
                    metadata.setdefault(fields[0], []).append([_clean(value) for value in fields[1:]])
                continue
            if line == "!series_matrix_table_end":
                break
            fields = line.split("\t")
            if not fields:
                continue
            probe = _clean(fields[0])
            symbols = probe_symbols.get(probe, set()).intersection(allowed_genes)
            if not symbols:
                continue
            values = np.array([_numeric(value) for value in fields[1 : len(sample_ids) + 1]], dtype=np.float32)
            if values.size != len(sample_ids):
                continue
            valid = np.isfinite(values)
            for gene in symbols:
                if gene not in gene_sum:
                    gene_sum[gene] = np.zeros(len(sample_ids), dtype=np.float64)
                    gene_count[gene] = np.zeros(len(sample_ids), dtype=np.float64)
                gene_sum[gene][valid] += values[valid]
                gene_count[gene][valid] += 1.0

    titles = metadata.get("!Sample_title", [[]])[0]
    chars = metadata.get("!Sample_characteristics_ch1", [])
    labels, ages, sexes = _parse_characteristics(chars, len(sample_ids))
    for idx, title in enumerate(titles[: len(sample_ids)]):
        lower = title.casefold()
        if re.search(r",\s*mci\s*,", lower):
            labels[idx] = "MCI"
        elif re.search(r",\s*ad\s*,", lower):
            labels[idx] = "AD"
        elif "control" in lower:
            labels[idx] = "HC"
    if include_unlabelled_rows:
        keep = list(range(len(sample_ids)))
    else:
        keep = [idx for idx, label in enumerate(labels) if label in {"HC", "MCI", "AD"}]
    expression = pd.DataFrame(
        {
            gene: np.divide(
                gene_sum[gene],
                gene_count[gene],
                out=np.full(len(sample_ids), np.nan, dtype=np.float64),
                where=gene_count[gene] > 0,
            )
            for gene in sorted(gene_sum)
        },
        index=sample_ids,
    )
    expression = expression.iloc[keep].astype(np.float32)
    meta = pd.DataFrame(
        {
            "label": [labels[idx] for idx in keep],
            "age": ages[keep],
            "sex_num": sexes[keep],
            "dataset": "GSE140829",
            "study_family": "GSE140829",
        },
        index=expression.index,
    )
    audit = {
        "series_sample_count": len(sample_ids),
        "retained_sample_count": len(keep),
        "label_counts": dict(sorted(meta["label"].value_counts().to_dict().items())),
        "platform_probe_count": len(probe_symbols),
        "mapped_gene_count": len(gene_sum),
        "missing_locked_symbols": [gene for gene, _ in PANEL if gene not in gene_sum],
        "locked_symbols_present": [gene for gene, _ in PANEL if gene in gene_sum],
    }
    return expression, meta, audit


def make_estimator(kind: str, feature_count: int, seed: int) -> Pipeline:
    if kind == "transcriptome":
        pipeline = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                ("select", SelectKBest(score_func=f_classif, k=min(200, feature_count))),
                ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
            ]
        )
    else:
        pipeline = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
            ]
        )
    return pipeline


def metric_row(y: np.ndarray, probability: np.ndarray, status: str, endpoint: str, cohort: str, model: str) -> dict[str, Any]:
    auc = float(roc_auc_score(y, probability)) if len(np.unique(y)) == 2 else np.nan
    brier = float(brier_score_loss(y, probability))
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    logit = np.log(clipped / (1 - clipped))
    try:
        calibration = LogisticRegression(C=1e6, solver="lbfgs", max_iter=1000)
        calibration.fit(logit.reshape(-1, 1), y)
        cal_intercept = float(calibration.intercept_[0])
        cal_slope = float(calibration.coef_[0, 0])
    except (ValueError, np.linalg.LinAlgError):
        cal_intercept, cal_slope = np.nan, np.nan
    bins = np.linspace(0, 1, 11)
    ece = 0.0
    used_bins = 0
    for low, high in zip(bins[:-1], bins[1:]):
        mask = (probability >= low) & (probability < high if high < 1 else probability <= high)
        if not mask.any():
            continue
        ece += float(mask.mean()) * abs(float(probability[mask].mean()) - float(y[mask].mean()))
        used_bins += 1
    return {
        "status": status,
        "endpoint": endpoint,
        "cohort": cohort,
        "model": model,
        "n": int(len(y)),
        "positive_n": int(y.sum()),
        "auc": auc,
        "brier": brier,
        "calibration_intercept": cal_intercept,
        "calibration_slope": cal_slope,
        "calibration_method": "descriptive_logistic_recalibration",
        "calibration_fit_scope": "evaluation_labels_descriptive_only",
        "ece": float(ece),
        "ece_nonempty_bins": used_bins,
    }


def bootstrap_auc(y: np.ndarray, probability: np.ndarray, seed: int, repetitions: int = 1000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    for _ in range(repetitions):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        estimates.append(float(roc_auc_score(y[idx], probability[idx])))
    if not estimates:
        return np.nan, np.nan
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def frozen_score(expression: pd.DataFrame) -> tuple[pd.Series, list[str]]:
    present = [gene for gene, _ in PANEL if gene in expression.columns]
    missing = [gene for gene, _ in PANEL if gene not in expression.columns]
    score = pd.Series(0.0, index=expression.index)
    for gene, sign in PANEL:
        if gene in expression.columns:
            score = score + sign * expression[gene].astype(float) / 6.0
    return score, missing


def run_endpoint(
    endpoint: str,
    dev_x: pd.DataFrame,
    dev_meta: pd.DataFrame,
    ext_x: pd.DataFrame,
    ext_meta: pd.DataFrame,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    if endpoint == "stage_mci_vs_hc":
        labels = {"MCI": 1, "HC": 0}
    elif endpoint == "stage_mci_vs_ad":
        labels = {"MCI": 1, "AD": 0}
    else:
        labels = {"AD": 1, "HC": 0}
    dev_keep = dev_meta["label"].isin(labels)
    ext_keep = ext_meta["label"].isin(labels)
    y_dev = dev_meta.loc[dev_keep, "label"].map(labels).to_numpy(dtype=int)
    y_ext = ext_meta.loc[ext_keep, "label"].map(labels).to_numpy(dtype=int)
    x_dev = dev_x.loc[dev_keep]
    x_ext = ext_x.loc[ext_keep]
    dev_cov = dev_meta.loc[dev_keep, ["age", "sex_num"]].to_numpy(dtype=float)
    ext_cov = ext_meta.loc[ext_keep, ["age", "sex_num"]].to_numpy(dtype=float)
    metrics: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    selections: list[dict[str, Any]] = []
    selective_curves: list[dict[str, Any]] = []
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    model_inputs = {
        "age_only": (dev_cov[:, [0]], ext_cov[:, [0]], 1),
        "sex_only": (dev_cov[:, [1]], ext_cov[:, [1]], 1),
        "age_plus_sex": (dev_cov, ext_cov, 2),
        "transcriptome": (x_dev, x_ext, x_dev.shape[1]),
    }
    for model_name, (train_input, external_input, feature_count) in model_inputs.items():
        oof = np.full(len(y_dev), np.nan, dtype=float)
        for fold, (train_idx, test_idx) in enumerate(outer.split(np.zeros(len(y_dev)), y_dev), start=1):
            estimator = make_estimator(model_name, feature_count, seed + fold)
            train_slice = train_input.iloc[train_idx] if hasattr(train_input, "iloc") else train_input[train_idx]
            test_slice = train_input.iloc[test_idx] if hasattr(train_input, "iloc") else train_input[test_idx]
            estimator.fit(train_slice, y_dev[train_idx])
            oof[test_idx] = estimator.predict_proba(test_slice)[:, 1]
            if model_name == "transcriptome":
                selected = estimator.named_steps["select"].get_support()
                for gene, keep in zip(train_input.columns, selected):
                    if keep:
                        selections.append({"endpoint": endpoint, "model": model_name, "fold": fold, "gene": gene})
        row = metric_row(y_dev, oof, "NESTED_INTERNAL", endpoint, "GSE63060_GSE63061", model_name)
        row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y_dev, oof, seed + 100)
        metrics.append(row)
        final_estimator = make_estimator(model_name, feature_count, seed + 999)
        final_estimator.fit(train_input, y_dev)
        probability = final_estimator.predict_proba(external_input)[:, 1]
        ext_row = metric_row(y_ext, probability, "EXTERNAL_VALIDATION", endpoint, "GSE140829", model_name)
        ext_row["auc_ci_low"], ext_row["auc_ci_high"] = bootstrap_auc(y_ext, probability, seed + 200)
        if model_name == "transcriptome":
            best = final_estimator
            train_transformed = best.named_steps["imputer"].transform(train_input)
            train_transformed = best.named_steps["scale"].transform(train_transformed)
            train_transformed = best.named_steps["select"].transform(train_transformed)
            ext_transformed = best.named_steps["imputer"].transform(external_input)
            ext_transformed = best.named_steps["scale"].transform(ext_transformed)
            ext_transformed = best.named_steps["select"].transform(ext_transformed)
            train_distance = np.nanmean(np.abs(train_transformed), axis=1)
            ext_distance = np.nanmean(np.abs(ext_transformed), axis=1)
            ood_threshold = float(np.quantile(train_distance, 0.95))
            ext_row["ood_rate"] = float(np.mean(ext_distance > ood_threshold))
            ext_row["ood_threshold_train_q95"] = ood_threshold
            ext_row["selected_gene_count"] = int(train_transformed.shape[1])
        else:
            ext_distance = np.full(len(y_ext), np.nan)
            ext_row["ood_rate"] = np.nan
            ext_row["ood_threshold_train_q95"] = np.nan
            ext_row["selected_gene_count"] = np.nan
        dev_confidence = np.maximum(oof, 1 - oof)
        confidence_threshold = float(np.quantile(dev_confidence, 0.20))
        ext_confidence = np.maximum(probability, 1 - probability)
        accepted = ext_confidence >= confidence_threshold
        ext_row["selective_threshold_train_q20"] = confidence_threshold
        ext_row["acceptance_rate"] = float(np.mean(accepted))
        ext_row["selective_risk"] = float(np.mean((probability[accepted] >= 0.5) != y_ext[accepted])) if accepted.any() else np.nan
        metrics.append(ext_row)
        for idx, (truth, pred, accept, distance) in enumerate(zip(y_ext, probability, accepted, ext_distance)):
            predictions.append(
                {
                    "endpoint": endpoint,
                    "cohort": "GSE140829",
                    "model": model_name,
                    "assessment_index": idx,
                    "y": int(truth),
                    "probability": float(pred),
                    "accepted": int(accept),
                    "ood_score": float(distance) if np.isfinite(distance) else np.nan,
                }
            )
        if endpoint == "stage_mci_vs_hc" and model_name == "transcriptome":
            ext_confidence = np.maximum(probability, 1 - probability)
            for quantile in np.linspace(0.0, 0.9, 10):
                development_threshold = float(np.quantile(dev_confidence, quantile))
                keep = ext_confidence >= development_threshold
                selective_curves.append(
                    {
                        "endpoint": endpoint,
                        "cohort": "GSE140829",
                        "model": model_name,
                        "development_confidence_quantile": float(quantile),
                        "development_threshold": development_threshold,
                        "external_coverage": float(np.mean(keep)),
                        "external_selective_error": float(
                            np.mean((probability[keep] >= 0.5) != y_ext[keep])
                        ) if keep.any() else np.nan,
                        "external_n": int(len(y_ext)),
                        "accepted_n": int(keep.sum()),
                    }
                )
    return metrics, predictions, selections, selective_curves


def add_frozen_score_rows(
    metrics: list[dict[str, Any]],
    dev_x: pd.DataFrame,
    dev_meta: pd.DataFrame,
    ext_x: pd.DataFrame,
    ext_meta: pd.DataFrame,
    seed: int,
) -> None:
    development_score, development_missing = frozen_score(dev_x)
    development_center = float(np.nanmedian(development_score))
    development_scale = float(np.nanstd(development_score) + 1e-6)
    for endpoint, labels in {
        "stage_mci_vs_hc": {"MCI": 1, "HC": 0},
        "stage_mci_vs_ad": {"MCI": 1, "AD": 0},
    }.items():
        for cohort, expression, meta in (
            ("GSE63060_GSE63061", dev_x, dev_meta),
            ("GSE140829", ext_x, ext_meta),
        ):
            score, missing = frozen_score(expression)
            keep = meta["label"].isin(labels)
            y = meta.loc[keep, "label"].map(labels).to_numpy(dtype=int)
            values = score.loc[keep].to_numpy(dtype=float)
            probabilities = 1.0 / (1.0 + np.exp(-((values - development_center) / development_scale)))
            row = metric_row(
                y,
                probabilities,
                "FROZEN_EXACT_INTERNAL" if cohort.startswith("GSE630") and not missing else "PARTIAL_EXTERNAL_SENSITIVITY",
                endpoint,
                cohort,
                "frozen_12_gene_score",
            )
            row["missing_locked_symbols"] = ";".join(missing)
            row["development_missing_locked_symbols"] = ";".join(development_missing)
            row["exact_mapping"] = int(len(missing) == 0)
            row["score_contrast_mean_positive_minus_negative"] = float(values[y == 1].mean() - values[y == 0].mean())
            row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y, values, seed + 300)
            metrics.append(row)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    dev_genes = set(dev60_x.columns).intersection(dev61_x.columns)
    dev_full_x = pd.concat([dev60_x[sorted(dev_genes)], dev61_x[sorted(dev_genes)]], axis=0)
    dev_meta = pd.concat([dev60_meta, dev61_meta], axis=0).loc[dev_full_x.index]
    ext_full_x, ext_meta, ext_audit = load_external(GSE140829_MATRIX, GPL15988_FULL, dev_genes)
    common = sorted(set(dev_full_x.columns).intersection(ext_full_x.columns))
    dev_x = dev_full_x[common]
    ext_x = ext_full_x[common]
    metrics: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    selections: list[dict[str, Any]] = []
    selective_curves: list[dict[str, Any]] = []
    for offset, endpoint in enumerate(["stage_mci_vs_hc", "stage_mci_vs_ad", "stage_ad_vs_hc_context"]):
        current_metrics, current_predictions, current_selections, current_curves = run_endpoint(
            endpoint, dev_x, dev_meta, ext_x, ext_meta, 20260812 + offset * 100
        )
        metrics.extend(current_metrics)
        predictions.extend(current_predictions)
        selections.extend(current_selections)
        selective_curves.extend(current_curves)
    add_frozen_score_rows(metrics, dev_full_x, dev_meta, ext_full_x, ext_meta, 20260812)
    pd.DataFrame(metrics).to_csv(OUT / "metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(OUT / "predictions_anonymized.csv", index=False)
    pd.DataFrame(selections).to_csv(OUT / "selection_events.csv", index=False)
    pd.DataFrame(selective_curves).to_csv(OUT / "selective_risk_curve.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "protocol_sha256": sha256_file(PROTOCOL),
        "script_sha256": sha256_file(Path(__file__)),
        "development": {
            "datasets": ["GSE63060", "GSE63061"],
            "sample_count": int(len(dev_meta)),
            "label_counts": dict(sorted(dev_meta["label"].value_counts().to_dict().items())),
            "gene_count_before_external_intersection": int(len(dev_genes)),
        },
        "external": ext_audit,
        "common_gene_count": int(len(common)),
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS, GSE140829_MATRIX, GPL15988_FULL]
        },
        "endpoints": ["stage_mci_vs_hc", "stage_mci_vs_ad", "stage_ad_vs_hc_context"],
        "external_label_access_used_only_for_assessment": True,
        "target_adaptation": False,
    }
    (OUT / "run_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"metrics": len(metrics), "predictions": len(predictions), "common_genes": len(common), "external": ext_audit}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
