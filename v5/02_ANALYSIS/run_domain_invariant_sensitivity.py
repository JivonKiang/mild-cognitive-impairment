from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline

from run_pairwise_external import (
    GSE140829_MATRIX,
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    GPL15988_FULL,
    bootstrap_auc,
    load_dev,
    load_external,
    metric_row,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260812


def sample_rank_transform(frame: pd.DataFrame) -> pd.DataFrame:
    """Transform each sample independently; no cohort-level target information."""
    values = frame.to_numpy(dtype=float)
    ranks = np.apply_along_axis(lambda row: pd.Series(row).rank(method="average").to_numpy(), 1, values)
    scale = np.nanstd(ranks, axis=1, keepdims=True)
    center = np.nanmean(ranks, axis=1, keepdims=True)
    transformed = (ranks - center) / np.where(scale > 0, scale, 1.0)
    return pd.DataFrame(transformed, index=frame.index, columns=frame.columns).astype(np.float32)


def estimator(feature_count: int, seed: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("select", SelectKBest(score_func=f_classif, k=min(200, feature_count))),
            ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
        ]
    )


def main() -> int:
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    dev_genes = sorted(set(dev60_x.columns).intersection(dev61_x.columns))
    dev_x = pd.concat([dev60_x[dev_genes], dev61_x[dev_genes]], axis=0)
    dev_meta = pd.concat([dev60_meta, dev61_meta], axis=0).loc[dev_x.index]
    ext_x, ext_meta, ext_audit = load_external(GSE140829_MATRIX, GPL15988_FULL, set(dev_genes))
    common = sorted(set(dev_x.columns).intersection(ext_x.columns))
    dev_x = dev_x[common]
    ext_x = ext_x[common]
    labels = {"HC": 0, "MCI": 1}
    dev_keep = dev_meta["label"].isin(labels)
    ext_keep = ext_meta["label"].isin(labels)
    x_dev = sample_rank_transform(dev_x.loc[dev_keep])
    x_ext = sample_rank_transform(ext_x.loc[ext_keep])
    y_dev = dev_meta.loc[dev_keep, "label"].map(labels).to_numpy(dtype=int)
    y_ext = ext_meta.loc[ext_keep, "label"].map(labels).to_numpy(dtype=int)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = np.full(len(y_dev), np.nan, dtype=float)
    for fold, (train_idx, test_idx) in enumerate(cv.split(x_dev, y_dev), start=1):
        model = estimator(x_dev.shape[1], SEED + fold)
        model.fit(x_dev.iloc[train_idx], y_dev[train_idx])
        oof[test_idx] = model.predict_proba(x_dev.iloc[test_idx])[:, 1]
    final = estimator(x_dev.shape[1], SEED + 999)
    final.fit(x_dev, y_dev)
    external_probability = final.predict_proba(x_ext)[:, 1]
    rows = []
    internal = metric_row(y_dev, oof, "NESTED_INTERNAL_SENSITIVITY", "stage_mci_vs_hc", "GSE63060_GSE63061", "sample_rank_transform")
    internal["auc_ci_low"], internal["auc_ci_high"] = bootstrap_auc(y_dev, oof, SEED + 100)
    internal["transformation"] = "within-sample rank then within-sample z-score"
    internal["target_labels_used_for_transform"] = False
    rows.append(internal)
    external = metric_row(y_ext, external_probability, "EXTERNAL_SENSITIVITY", "stage_mci_vs_hc", "GSE140829", "sample_rank_transform")
    external["auc_ci_low"], external["auc_ci_high"] = bootstrap_auc(y_ext, external_probability, SEED + 200)
    external["transformation"] = "within-sample rank then within-sample z-score"
    external["target_labels_used_for_transform"] = False
    external["source_external_label_access"] = "assessment only"
    rows.append(external)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "domain_invariant_sensitivity.csv", index=False)
    pd.DataFrame(
        {
            "cohort": "GSE140829",
            "assessment_index": np.arange(len(y_ext)),
            "y": y_ext,
            "probability": external_probability,
        }
    ).to_csv(OUT / "domain_invariant_predictions.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "design": "sample-wise rank/z-score sensitivity; no cohort-level target adaptation",
        "external_labels_used_only_for_assessment": True,
        "external_audit": ext_audit,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS, GSE140829_MATRIX, GPL15988_FULL]
        },
    }
    (OUT / "domain_invariant_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(pd.DataFrame(rows).to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
