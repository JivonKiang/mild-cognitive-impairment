from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit

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
    make_estimator,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260812
REPEATS = 20
ALPHA = 0.10


def conformal_set(probability: np.ndarray, q: float) -> list[str]:
    sets = []
    for p in probability:
        members = []
        if 1.0 - (1.0 - p) <= q:
            members.append("HC")
        if 1.0 - p <= q:
            members.append("MCI")
        sets.append("|".join(members) if members else "empty")
    return sets


def quantile_higher(values: np.ndarray, alpha: float) -> float:
    n = len(values)
    rank = int(np.ceil((n + 1) * (1 - alpha))) - 1
    rank = min(max(rank, 0), n - 1)
    return float(np.sort(values)[rank])


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
    x_dev, y_dev = x_dev.loc[dev_keep], m_dev.loc[dev_keep, "label"].map(labels).to_numpy(dtype=int)
    x_ext, y_ext = x_ext.loc[ext_keep], m_ext.loc[ext_keep, "label"].map(labels).to_numpy(dtype=int)
    rows = []
    predictions = []
    splitter = StratifiedShuffleSplit(n_splits=REPEATS, test_size=0.30, random_state=SEED)
    for repeat, (proper_idx, calibration_idx) in enumerate(splitter.split(x_dev, y_dev), start=1):
        model = make_estimator("transcriptome", x_dev.shape[1], SEED + repeat)
        model.fit(x_dev.iloc[proper_idx], y_dev[proper_idx])
        calibration_probability = model.predict_proba(x_dev.iloc[calibration_idx])[:, 1]
        calibration_nonconformity = 1.0 - np.where(y_dev[calibration_idx] == 1, calibration_probability, 1.0 - calibration_probability)
        q = quantile_higher(calibration_nonconformity, ALPHA)
        external_probability = model.predict_proba(x_ext)[:, 1]
        sets = conformal_set(external_probability, q)
        included = [
            ("HC" in prediction_set and truth == 0) or ("MCI" in prediction_set and truth == 1)
            for prediction_set, truth in zip(sets, y_ext)
        ]
        singleton = np.array([prediction_set in {"HC", "MCI"} for prediction_set in sets])
        singleton_error = np.mean(
            [
                ((prediction_set == "MCI") != bool(truth))
                for prediction_set, truth in zip(sets, y_ext)
                if prediction_set in {"HC", "MCI"}
            ]
        ) if singleton.any() else np.nan
        rows.append(
            {
                "repeat": repeat,
                "alpha": ALPHA,
                "calibration_n": int(len(calibration_idx)),
                "proper_training_n": int(len(proper_idx)),
                "external_n": int(len(y_ext)),
                "q": q,
                "external_coverage": float(np.mean(included)),
                "mean_prediction_set_size": float(np.mean([0 if s == "empty" else len(s.split("|")) for s in sets])),
                "ambiguous_rate": float(np.mean([s == "HC|MCI" for s in sets])),
                "empty_rate": float(np.mean([s == "empty" for s in sets])),
                "singleton_rate": float(np.mean(singleton)),
                "singleton_error": float(singleton_error) if np.isfinite(singleton_error) else np.nan,
                "external_auc": float(__import__("sklearn.metrics").metrics.roc_auc_score(y_ext, external_probability)),
            }
        )
        if repeat == 1:
            for index, (truth, probability, prediction_set, covered) in enumerate(zip(y_ext, external_probability, sets, included)):
                predictions.append(
                    {
                        "repeat": repeat,
                        "assessment_index": index,
                        "y": int(truth),
                        "probability": float(probability),
                        "prediction_set": prediction_set,
                        "covered": int(covered),
                    }
                )
    summary = pd.DataFrame(rows)
    aggregate = pd.DataFrame(
        [
            {
                "alpha": ALPHA,
                "repeats": REPEATS,
                "external_coverage_mean": summary.external_coverage.mean(),
                "external_coverage_q025": summary.external_coverage.quantile(0.025),
                "external_coverage_q975": summary.external_coverage.quantile(0.975),
                "mean_prediction_set_size": summary.mean_prediction_set_size.mean(),
                "ambiguous_rate_mean": summary.ambiguous_rate.mean(),
                "singleton_rate_mean": summary.singleton_rate.mean(),
                "singleton_error_mean": summary.singleton_error.mean(),
                "external_auc_mean": summary.external_auc.mean(),
                "external_auc_q025": summary.external_auc.quantile(0.025),
                "external_auc_q975": summary.external_auc.quantile(0.975),
            }
        ]
    )
    OUT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT / "split_conformal_replicates.csv", index=False)
    aggregate.to_csv(OUT / "split_conformal_summary.csv", index=False)
    pd.DataFrame(predictions).to_csv(OUT / "split_conformal_predictions_repeat1.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "design": "split conformal prediction using development-only calibration; external labels assessment-only",
        "alpha": ALPHA,
        "repeats": REPEATS,
        "external_labels_used_for_fit": False,
        "target_domain_calibration_or_adaptation": False,
        "external_audit": ext_audit,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS, GSE140829_MATRIX, GPL15988_FULL]
        },
    }
    (OUT / "split_conformal_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(aggregate.to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
