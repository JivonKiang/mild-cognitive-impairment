from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from run_pairwise_external import (
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    PANEL,
    bootstrap_auc,
    frozen_score,
    load_dev,
    make_estimator,
    metric_row,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"


def endpoint_view(
    expression: pd.DataFrame, metadata: pd.DataFrame, labels: dict[str, int]
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    keep = metadata["label"].isin(labels)
    x = expression.loc[keep]
    m = metadata.loc[keep]
    y = m["label"].map(labels).to_numpy(dtype=int)
    return x, m, y


def transformed_for_ood(estimator: Any, values: pd.DataFrame) -> np.ndarray:
    transformed = estimator.named_steps["imputer"].transform(values)
    transformed = estimator.named_steps["scale"].transform(transformed)
    if "select" in estimator.named_steps:
        transformed = estimator.named_steps["select"].transform(transformed)
    return transformed


def run_direction(
    train_id: str,
    test_id: str,
    train_x: pd.DataFrame,
    train_meta: pd.DataFrame,
    test_x: pd.DataFrame,
    test_meta: pd.DataFrame,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    labels = {"MCI": 1, "HC": 0}
    train_x, train_meta, y_train = endpoint_view(train_x, train_meta, labels)
    test_x, test_meta, y_test = endpoint_view(test_x, test_meta, labels)
    train_cov = train_meta[["age", "sex_num"]].to_numpy(dtype=float)
    test_cov = test_meta[["age", "sex_num"]].to_numpy(dtype=float)
    model_inputs = {
        "age_only": (train_cov[:, [0]], test_cov[:, [0]], 1),
        "sex_only": (train_cov[:, [1]], test_cov[:, [1]], 1),
        "age_plus_sex": (train_cov, test_cov, 2),
        "transcriptome": (train_x, test_x, train_x.shape[1]),
    }
    metrics: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    selections: list[dict[str, Any]] = []
    for offset, (model_name, (x_train, x_test, feature_count)) in enumerate(model_inputs.items()):
        estimator = make_estimator(model_name, feature_count, seed + offset)
        estimator.fit(x_train, y_train)
        probability = estimator.predict_proba(x_test)[:, 1]
        row = metric_row(
            y_test,
            probability,
            "STUDY_HELD_OUT",
            "stage_mci_vs_hc",
            f"{train_id}_TO_{test_id}",
            model_name,
        )
        row["train_domain"] = train_id
        row["test_domain"] = test_id
        row["train_n"] = int(len(y_train))
        row["test_n"] = int(len(y_test))
        row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y_test, probability, seed + 100)
        row["target_labels_used_for_fit"] = False
        if model_name == "transcriptome":
            transformed_train = transformed_for_ood(estimator, x_train)
            transformed_test = transformed_for_ood(estimator, x_test)
            train_distance = np.nanmean(np.abs(transformed_train), axis=1)
            test_distance = np.nanmean(np.abs(transformed_test), axis=1)
            threshold = float(np.quantile(train_distance, 0.95))
            row["ood_rate"] = float(np.mean(test_distance > threshold))
            row["ood_threshold_train_q95"] = threshold
            row["selected_gene_count"] = int(transformed_train.shape[1])
            support = estimator.named_steps["select"].get_support()
            selected_genes = [gene for gene, selected in zip(x_train.columns, support) if selected]
            for gene in selected_genes:
                selections.append(
                    {
                        "train_domain": train_id,
                        "test_domain": test_id,
                        "endpoint": "stage_mci_vs_hc",
                        "model": model_name,
                        "gene": gene,
                    }
                )
        else:
            test_distance = np.full(len(y_test), np.nan)
            row["ood_rate"] = np.nan
            row["ood_threshold_train_q95"] = np.nan
            row["selected_gene_count"] = np.nan
        metrics.append(row)
        for index, (truth, pred, distance) in enumerate(zip(y_test, probability, test_distance)):
            predictions.append(
                {
                    "endpoint": "stage_mci_vs_hc",
                    "train_domain": train_id,
                    "test_domain": test_id,
                    "model": model_name,
                    "assessment_index": index,
                    "y": int(truth),
                    "probability": float(pred),
                    "ood_score": float(distance) if np.isfinite(distance) else np.nan,
                }
            )
    return metrics, predictions, selections


def frozen_transport_row(
    train_id: str,
    test_id: str,
    train_x: pd.DataFrame,
    train_meta: pd.DataFrame,
    test_x: pd.DataFrame,
    test_meta: pd.DataFrame,
    seed: int,
) -> dict[str, Any]:
    labels = {"MCI": 1, "HC": 0}
    train_keep = train_meta["label"].isin(labels)
    test_keep = test_meta["label"].isin(labels)
    y = test_meta.loc[test_keep, "label"].map(labels).to_numpy(dtype=int)
    train_x = train_x.loc[train_keep]
    test_x = test_x.loc[test_keep]
    present = [gene for gene, _sign in PANEL if gene in train_x.columns and gene in test_x.columns]
    missing_train = [gene for gene, _sign in PANEL if gene not in train_x.columns]
    missing_test = [gene for gene, _sign in PANEL if gene not in test_x.columns]
    train_standardized = pd.DataFrame(index=train_x.index)
    test_standardized = pd.DataFrame(index=test_x.index)
    for gene in present:
        center = float(train_x[gene].median())
        scale = float(train_x[gene].std(ddof=0) + 1e-6)
        train_standardized[gene] = (train_x[gene] - center) / scale
        test_standardized[gene] = (test_x[gene] - center) / scale
    train_score = pd.Series(0.0, index=train_x.index)
    test_score = pd.Series(0.0, index=test_x.index)
    for gene, sign in PANEL:
        if gene in present:
            train_score = train_score + sign * train_standardized[gene]
            test_score = test_score + sign * test_standardized[gene]
    values = test_score.to_numpy(dtype=float)
    auc = float(roc_auc_score(y, values))
    low, high = bootstrap_auc(y, values, seed + 300)
    return {
        "status": "FROZEN_STUDY_HELD_OUT",
        "endpoint": "stage_mci_vs_hc",
        "cohort": f"{train_id}_TO_{test_id}",
        "model": "frozen_12_gene_rank_score",
        "train_domain": train_id,
        "test_domain": test_id,
        "n": int(len(y)),
        "positive_n": int(y.sum()),
        "auc": auc,
        "auc_ci_low": low,
        "auc_ci_high": high,
        "brier": np.nan,
        "calibration_intercept": np.nan,
        "calibration_slope": np.nan,
        "ece": np.nan,
        "ood_rate": np.nan,
        "missing_locked_symbols_train": ";".join(missing_train),
        "missing_locked_symbols_test": ";".join(missing_test),
        "present_locked_symbol_count": len(present),
        "exact_mapping": int(not missing_train and not missing_test),
        "score_contrast_mean_positive_minus_negative": float(
            test_score[y == 1].mean() - test_score[y == 0].mean()
        ),
        "target_labels_used_for_fit": False,
    }


def main() -> int:
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    common = sorted(set(dev60_x.columns).intersection(dev61_x.columns))
    dev60_x = dev60_x[common]
    dev61_x = dev61_x[common]
    metrics: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    selections: list[dict[str, Any]] = []
    for offset, (train_id, test_id, train_x, train_meta, test_x, test_meta) in enumerate(
        [
            ("GSE63060", "GSE63061", dev60_x, dev60_meta, dev61_x, dev61_meta),
            ("GSE63061", "GSE63060", dev61_x, dev61_meta, dev60_x, dev60_meta),
        ]
    ):
        current_metrics, current_predictions, current_selections = run_direction(
            train_id, test_id, train_x, train_meta, test_x, test_meta, 20260812 + offset * 1000
        )
        metrics.extend(current_metrics)
        predictions.extend(current_predictions)
        selections.extend(current_selections)
        metrics.append(frozen_transport_row(train_id, test_id, train_x, train_meta, test_x, test_meta, 20260812 + offset))
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metrics).to_csv(OUT / "study_heldout_metrics.csv", index=False)
    pd.DataFrame(predictions).to_csv(OUT / "study_heldout_predictions.csv", index=False)
    pd.DataFrame(selections).to_csv(OUT / "study_heldout_selection_events.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "endpoint": "stage_mci_vs_hc",
        "design": "bidirectional study-held-out transportability",
        "domains": ["GSE63060", "GSE63061"],
        "common_exact_symbol_gene_count": len(common),
        "feature_selection_inside_training_domain_only": True,
        "external_test_labels_used_for_fit": False,
        "target_domain_calibration_or_adaptation": False,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS]
        },
    }
    (OUT / "study_heldout_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"metrics": len(metrics), "predictions": len(predictions), "common_genes": len(common)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
