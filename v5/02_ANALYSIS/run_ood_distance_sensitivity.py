from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

from run_pairwise_external import (
    GSE140829_MATRIX,
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    GPL15988_FULL,
    load_dev,
    load_external,
    make_estimator,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
PROTOCOL = ROOT / "00_PROTOCOL" / "analysis_protocol_pairwise_v1.md"
SEED = 20260812
THRESHOLDS = (0.90, 0.95, 0.99)


def transform_locked_pipeline(model, frame: pd.DataFrame) -> np.ndarray:
    """Apply the already-fitted development pipeline without refitting."""
    values = model.named_steps["imputer"].transform(frame)
    values = model.named_steps["scale"].transform(values)
    return model.named_steps["select"].transform(values)


def mean_abs_distance(train_values: np.ndarray, test_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.nanmean(np.abs(train_values), axis=1), np.nanmean(np.abs(test_values), axis=1)


def mahalanobis_ledoitwolf(train_values: np.ndarray, test_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    estimator = LedoitWolf().fit(train_values)
    train_centered = train_values - estimator.location_
    test_centered = test_values - estimator.location_
    train_distance = np.sqrt(np.maximum(0.0, np.einsum("ij,jk,ik->i", train_centered, estimator.precision_, train_centered)))
    test_distance = np.sqrt(np.maximum(0.0, np.einsum("ij,jk,ik->i", test_centered, estimator.precision_, test_centered)))
    return train_distance, test_distance


def robust_pca_distance(train_values: np.ndarray, test_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Label-free PCA distance using training-only median/MAD score scaling."""
    n_components = min(20, train_values.shape[0] - 1, train_values.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full", random_state=SEED)
    train_scores = pca.fit_transform(train_values)
    test_scores = pca.transform(test_values)
    center = np.median(train_scores, axis=0)
    mad = np.median(np.abs(train_scores - center), axis=0) * 1.4826
    scale = np.where(mad > 1e-12, mad, 1.0)
    train_scaled = (train_scores - center) / scale
    test_scaled = (test_scores - center) / scale
    return np.sqrt(np.mean(train_scaled**2, axis=1)), np.sqrt(np.mean(test_scaled**2, axis=1))


def rate(values: np.ndarray) -> float:
    return float(np.mean(values)) if len(values) else float("nan")


def main() -> int:
    x60, m60 = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    x61, m61 = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    genes = sorted(set(x60.columns).intersection(x61.columns))
    x_dev = pd.concat([x60[genes], x61[genes]], axis=0)
    m_dev = pd.concat([m60, m61], axis=0).loc[x_dev.index]
    x_ext, m_ext, ext_audit = load_external(GSE140829_MATRIX, GPL15988_FULL, set(genes))
    common = sorted(set(x_dev.columns).intersection(x_ext.columns))
    x_dev, x_ext = x_dev[common], x_ext[common]
    labels = {"HC": 0, "MCI": 1}
    dev_keep = m_dev["label"].isin(labels)
    ext_keep = m_ext["label"].isin(labels)
    x_dev = x_dev.loc[dev_keep]
    x_ext = x_ext.loc[ext_keep]
    y_dev = m_dev.loc[dev_keep, "label"].map(labels).to_numpy(dtype=int)
    y_ext = m_ext.loc[ext_keep, "label"].map(labels).to_numpy(dtype=int)

    model = make_estimator("transcriptome", x_dev.shape[1], SEED + 999)
    model.fit(x_dev, y_dev)
    train_values = transform_locked_pipeline(model, x_dev)
    ext_values = transform_locked_pipeline(model, x_ext)
    probabilities = model.predict_proba(x_ext)[:, 1]
    errors = (probabilities >= 0.5) != y_ext

    distance_functions = {
        "mean_abs_standardized": mean_abs_distance,
        "mahalanobis_ledoitwolf": mahalanobis_ledoitwolf,
        "robust_pca_distance": robust_pca_distance,
    }
    rows: list[dict[str, object]] = []
    for method, function in distance_functions.items():
        train_distance, ext_distance = function(train_values, ext_values)
        for quantile in THRESHOLDS:
            threshold = float(np.quantile(train_distance, quantile))
            ood = ext_distance > threshold
            inlier = ~ood
            rows.append(
                {
                    "endpoint": "stage_mci_vs_hc",
                    "cohort": "GSE140829",
                    "distance_method": method,
                    "threshold_quantile": quantile,
                    "threshold_train": threshold,
                    "development_n": int(len(y_dev)),
                    "external_n": int(len(y_ext)),
                    "external_ood_n": int(ood.sum()),
                    "external_inlier_n": int(inlier.sum()),
                    "external_ood_rate": rate(ood),
                    "external_error_rate_ood": rate(errors[ood]),
                    "external_error_rate_inlier": rate(errors[inlier]),
                    "external_error_rate_difference_ood_minus_inlier": (
                        rate(errors[ood]) - rate(errors[inlier])
                        if ood.any() and inlier.any()
                        else float("nan")
                    ),
                    "external_auc_assessment_only": float(roc_auc_score(y_ext, probabilities)),
                    "external_labels_used_for_fit": False,
                    "external_labels_used_for_threshold": False,
                }
            )
    result = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / "ood_distance_sensitivity.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "design": "training-only OOD distance and threshold sensitivity audit",
        "endpoint": "stage_mci_vs_hc",
        "distance_methods": sorted(distance_functions),
        "threshold_quantiles": list(THRESHOLDS),
        "external_labels_used_for_fit": False,
        "external_labels_used_for_threshold": False,
        "external_labels_used_only_for_error_and_auc_assessment": True,
        "target_domain_adaptation": False,
        "external_audit": ext_audit,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS, GSE140829_MATRIX, GPL15988_FULL]
        },
        "protocol_sha256": sha256_file(PROTOCOL),
    }
    (OUT / "ood_distance_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(result.to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
