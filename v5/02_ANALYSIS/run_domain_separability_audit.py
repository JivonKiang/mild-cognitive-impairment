from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from run_pairwise_external import (
    GSE140829_MATRIX,
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    GPL15988_FULL,
    load_dev,
    load_external,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260812
REPEATS = 10
FOLDS = 5
PERMUTATIONS = 5000


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def estimator(feature_count: int, seed: int) -> Pipeline:
    # Domain labels are intentionally the only labels used by this audit.
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("select", SelectKBest(score_func=f_classif, k=min(200, feature_count))),
            ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
        ]
    )


def bootstrap_metrics(y: np.ndarray, probability: np.ndarray, seed: int, repetitions: int = 2000) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    auc_values: list[float] = []
    balanced_values: list[float] = []
    for _ in range(repetitions):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        auc_values.append(float(roc_auc_score(y[idx], probability[idx])))
        balanced_values.append(float(balanced_accuracy_score(y[idx], probability[idx] >= 0.5)))
    return {
        "auc_ci_low": float(np.quantile(auc_values, 0.025)),
        "auc_ci_high": float(np.quantile(auc_values, 0.975)),
        "balanced_accuracy_ci_low": float(np.quantile(balanced_values, 0.025)),
        "balanced_accuracy_ci_high": float(np.quantile(balanced_values, 0.975)),
    }


def main() -> int:
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060", include_unlabelled_rows=True)
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061", include_unlabelled_rows=True)
    dev_genes = sorted(set(dev60_x.columns).intersection(dev61_x.columns))
    dev_x = pd.concat([dev60_x[dev_genes], dev61_x[dev_genes]], axis=0)
    ext_x, ext_meta, ext_audit = load_external(
        GSE140829_MATRIX, GPL15988_FULL, set(dev_genes), include_unlabelled_rows=True
    )
    common = sorted(set(dev_x.columns).intersection(ext_x.columns))
    source_x = dev_x[common].copy()
    external_x = ext_x[common].copy()
    expression = pd.concat([source_x, external_x], axis=0, ignore_index=True)
    domain_y = np.concatenate([
        np.zeros(len(source_x), dtype=int),
        np.ones(len(external_x), dtype=int),
    ])
    source_labels = pd.concat([dev60_meta, dev61_meta], axis=0).loc[source_x.index, "label"].to_numpy()
    external_labels = ext_meta.loc[external_x.index, "label"].to_numpy()

    # Disease labels are retained only for post hoc stratum sensitivities.
    # The all-row audit below does not use them for filtering, preprocessing,
    # feature selection, model fitting, thresholding or the permutation null.
    cv_predictions = np.full((REPEATS, len(domain_y)), np.nan, dtype=float)
    selection_counts: dict[str, int] = {}
    for repeat in range(REPEATS):
        cv = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=SEED + repeat)
        for fold, (train_idx, test_idx) in enumerate(cv.split(expression, domain_y), start=1):
            model = estimator(expression.shape[1], SEED + repeat * 100 + fold)
            model.fit(expression.iloc[train_idx], domain_y[train_idx])
            cv_predictions[repeat, test_idx] = model.predict_proba(expression.iloc[test_idx])[:, 1]
            selected = model.named_steps["select"].get_support()
            for gene, keep in zip(expression.columns, selected):
                if keep:
                    selection_counts[gene] = selection_counts.get(gene, 0) + 1

    probability = np.nanmean(cv_predictions, axis=0)
    auc = float(roc_auc_score(domain_y, probability))
    balanced_accuracy = float(balanced_accuracy_score(domain_y, probability >= 0.5))
    intervals = bootstrap_metrics(domain_y, probability, SEED + 1)
    replicate_rows = []
    for repeat in range(REPEATS):
        repeat_probability = cv_predictions[repeat]
        replicate_rows.append(
            {
                "repeat": repeat + 1,
                "auc": float(roc_auc_score(domain_y, repeat_probability)),
                "balanced_accuracy": float(
                    balanced_accuracy_score(domain_y, repeat_probability >= 0.5)
                ),
                "source_n": len(source_x),
                "external_n": len(external_x),
            }
        )

    rng = np.random.default_rng(SEED + 2)
    permutation_auc = np.empty(PERMUTATIONS, dtype=float)
    for idx in range(PERMUTATIONS):
        permutation_auc[idx] = roc_auc_score(rng.permutation(domain_y), probability)
    permutation_p = float((np.count_nonzero(permutation_auc >= auc) + 1) / (PERMUTATIONS + 1))

    predictions = pd.DataFrame(
        {
            "assessment_index": np.arange(len(domain_y)),
            "domain": np.where(domain_y == 1, "GSE140829_external", "AddNeuroMed_development"),
            "domain_label": domain_y,
            "domain_probability": probability,
                "source_rows": [int(len(source_x))] * len(domain_y),
            "external_rows": [int(len(external_x))] * len(domain_y),
        }
    )
    OUT.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(OUT / "domain_separability_predictions.csv", index=False)
    pd.DataFrame(replicate_rows).to_csv(OUT / "domain_separability_replicates.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "status": "DOMAIN_SEPARABILITY_AUDIT",
                "domain_comparison": "AddNeuroMed_development_vs_GSE140829_external",
                "source_n": len(source_x),
                "external_n": len(external_x),
                "common_gene_n": len(common),
                "auc": auc,
                "auc_ci_low": intervals["auc_ci_low"],
                "auc_ci_high": intervals["auc_ci_high"],
                "balanced_accuracy": balanced_accuracy,
                "balanced_accuracy_ci_low": intervals["balanced_accuracy_ci_low"],
                "balanced_accuracy_ci_high": intervals["balanced_accuracy_ci_high"],
                "permutation_n": PERMUTATIONS,
                "permutation_p": permutation_p,
                "cv_repeats": REPEATS,
                "cv_folds": FOLDS,
                "feature_selection_k": min(200, len(common)),
                "domain_labels_used_for_fit": True,
                "disease_labels_used_for_fit": False,
                "disease_labels_used_for_preprocessing": False,
                "disease_labels_used_for_row_alignment": True,
                "disease_label_filter_applied_by_loader": False,
                "disease_labels_used_for_stratum_definition": True,
                "external_labels_used_for_threshold": False,
                "analysis_unit": "sample row; GSE140829 subject independence unresolved",
            }
        ]
    )
    summary.to_csv(OUT / "domain_separability_summary.csv", index=False)
    top_features = pd.DataFrame(
        sorted(selection_counts.items(), key=lambda item: (-item[1], item[0])),
        columns=["gene", "selection_count_across_cv_fits"],
    )
    top_features.to_csv(OUT / "domain_separability_selected_features.csv", index=False)
    strata: list[dict[str, object]] = []
    for stratum in ["HC", "MCI", "AD"]:
        mask = np.concatenate([source_labels == stratum, external_labels == stratum])
        if mask.sum() < 20 or len(np.unique(domain_y[mask])) < 2:
            continue
        stratum_x = expression.loc[mask].reset_index(drop=True)
        stratum_y = domain_y[mask]
        stratum_probability = np.full(len(stratum_y), np.nan, dtype=float)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED + len(stratum))
        for fold, (train_idx, test_idx) in enumerate(cv.split(stratum_x, stratum_y), start=1):
            model = estimator(stratum_x.shape[1], SEED + len(stratum) * 100 + fold)
            model.fit(stratum_x.iloc[train_idx], stratum_y[train_idx])
            stratum_probability[test_idx] = model.predict_proba(stratum_x.iloc[test_idx])[:, 1]
        strata.append(
            {
                "stratum": stratum,
                "source_n": int(np.sum(source_labels == stratum)),
                "external_n": int(np.sum(external_labels == stratum)),
                "auc": float(roc_auc_score(stratum_y, stratum_probability)),
                "balanced_accuracy": float(balanced_accuracy_score(stratum_y, stratum_probability >= 0.5)),
                "analysis_role": "post hoc phenotype-stratified domain sensitivity",
                "disease_labels_used_for_domain_fit": False,
                "disease_labels_used_for_stratum_definition": True,
            }
        )
    pd.DataFrame(strata).to_csv(OUT / "domain_separability_strata.csv", index=False)
    audit = {
        "status": "RESULTS_GENERATED",
        "domain_label_definition": "0=AddNeuroMed development; 1=GSE140829 external",
        "target_disease_labels_used_for_fit": False,
        "target_disease_labels_used_for_preprocessing": False,
        "disease_labels_used_for_row_alignment": True,
        "disease_label_filter_applied_by_loader": False,
        "disease_labels_used_for_stratum_definition": True,
        "disease_label_filter_applied": False,
        "domain_labels_used_for_fit": True,
        "domain_labels_used_for_threshold": False,
        "feature_selection_scope": "within each training fold using domain labels only",
        "analysis_unit": "sample row; subject independence unresolved in GSE140829",
        "external_metadata_audit": ext_audit,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [GSE60_MATRIX, GSE60_LABELS, GSE61_MATRIX, GSE61_LABELS, GSE140829_MATRIX, GPL15988_FULL]
        },
        "script_sha256": sha256_file(Path(__file__)),
        "software": {"python": __import__("platform").python_version(), "scikit_learn": sklearn.__version__},
    }
    (OUT / "domain_separability_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(summary.to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
