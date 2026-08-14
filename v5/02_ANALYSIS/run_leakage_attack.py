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
from sklearn.preprocessing import StandardScaler

from run_pairwise_external import (
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    bootstrap_auc,
    load_dev,
    make_estimator,
    metric_row,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260812
REPEATS = 20


def load_primary() -> tuple[pd.DataFrame, pd.Series]:
    x60, m60 = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    x61, m61 = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    genes = sorted(set(x60.columns).intersection(x61.columns))
    x = pd.concat([x60[genes], x61[genes]], axis=0)
    meta = pd.concat([m60, m61], axis=0).loc[x.index]
    keep = meta["label"].isin({"HC", "MCI"})
    return x.loc[keep], meta.loc[keep, "label"].map({"HC": 0, "MCI": 1})


def attack_probability(x: pd.DataFrame, y: np.ndarray, seed: int) -> np.ndarray:
    # Deliberate attack: the label-aware selector sees every development sample
    # before cross-validation. This is a diagnostic comparator, not a valid model.
    preprocess = Pipeline(
        [("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler())]
    )
    x_all = preprocess.fit_transform(x)
    selector = SelectKBest(score_func=f_classif, k=min(200, x_all.shape[1]))
    x_selected = selector.fit_transform(x_all, y)
    oof = np.full(len(y), np.nan, dtype=float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for train_idx, test_idx in cv.split(x_selected, y):
        model = LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)
        model.fit(x_selected[train_idx], y[train_idx])
        oof[test_idx] = model.predict_proba(x_selected[test_idx])[:, 1]
    return oof


def clean_probability(x: pd.DataFrame, y: np.ndarray, seed: int) -> np.ndarray:
    oof = np.full(len(y), np.nan, dtype=float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (train_idx, test_idx) in enumerate(cv.split(x, y), start=1):
        model = make_estimator("transcriptome", x.shape[1], seed + fold)
        model.fit(x.iloc[train_idx], y[train_idx])
        oof[test_idx] = model.predict_proba(x.iloc[test_idx])[:, 1]
    return oof


def main() -> int:
    x, y_series = load_primary()
    y = y_series.to_numpy(dtype=int)
    rows: list[dict[str, object]] = []
    replicate_rows: list[dict[str, object]] = []
    for repeat in range(REPEATS):
        seed = SEED + repeat * 17
        for kind, probability in (
            ("nested_clean", clean_probability(x, y, seed)),
            ("full_scope_leakage_attack", attack_probability(x, y, seed)),
        ):
            row = metric_row(
                y,
                probability,
                "REPEATED_INTERNAL_DIAGNOSTIC",
                "stage_mci_vs_hc",
                "GSE63060_GSE63061",
                kind,
            )
            row["repeat"] = repeat + 1
            row["selection_scope"] = "training_fold" if kind == "nested_clean" else "full_development_before_cv"
            row["target_labels_used_for_external_fit"] = False
            row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y, probability, seed + 1000)
            rows.append(row)
            for idx, (truth, pred) in enumerate(zip(y, probability)):
                replicate_rows.append(
                    {
                        "repeat": repeat + 1,
                        "model": kind,
                        "assessment_index": idx,
                        "y": int(truth),
                        "probability": float(pred),
                    }
                )
    metrics = pd.DataFrame(rows)
    summary = (
        metrics.groupby("model", as_index=False)["auc"]
        .agg(["mean", "std", "min", "max"])
        .reset_index()
        .rename(columns={"mean": "auc_mean", "std": "auc_sd", "min": "auc_min", "max": "auc_max"})
    )
    summary["auc_q025"] = metrics.groupby("model")["auc"].quantile(0.025).to_numpy()
    summary["auc_q975"] = metrics.groupby("model")["auc"].quantile(0.975).to_numpy()
    summary["optimism_vs_nested_mean"] = np.nan
    clean_mean = float(summary.loc[summary["model"] == "nested_clean", "auc_mean"].iloc[0])
    summary.loc[summary["model"] == "full_scope_leakage_attack", "optimism_vs_nested_mean"] = (
        summary.loc[summary["model"] == "full_scope_leakage_attack", "auc_mean"] - clean_mean
    )
    OUT.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUT / "leakage_attack_replicates.csv", index=False)
    summary.to_csv(OUT / "leakage_attack_summary.csv", index=False)
    pd.DataFrame(replicate_rows).to_csv(OUT / "leakage_attack_predictions.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "endpoint": "stage_mci_vs_hc",
        "development_samples": int(len(y)),
        "label_counts": {"HC": int((y == 0).sum()), "MCI": int((y == 1).sum())},
        "repeats": REPEATS,
        "folds_per_repeat": 5,
        "clean_selection_scope": "training_fold_only",
        "attack_selection_scope": "full_development_before_cv",
        "attack_is_diagnostic_not_valid_model": True,
        "external_labels_used": False,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS]
        },
    }
    (OUT / "leakage_attack_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(summary.to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
