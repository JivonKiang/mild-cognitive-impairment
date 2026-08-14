from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from run_gse249477_rank_transport import (
    GSE249477_RAW,
    GSE249477_SERIES,
    load_target_expression,
    parse_assessment_labels as parse_gse249477_assessment_labels,
)
from run_pairwise_external import (
    GSE140829_MATRIX,
    GPL15988_FULL,
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    _clean,
    _numeric,
    bootstrap_auc,
    load_dev,
    load_external,
    load_gpl_symbol_map,
    metric_row,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260813
TOP_K = 100
OUTER_FOLDS = 5
PERMUTATION_REPETITIONS = 100
LABELS = {"HC": 0, "MCI": 1}


def study_effects(expression: pd.DataFrame, y: np.ndarray) -> pd.Series:
    """Return MCI-minus-HC standardized effects from supplied training rows only."""
    values = expression.to_numpy(dtype=float)
    cases = values[y == 1]
    controls = values[y == 0]
    if len(cases) < 2 or len(controls) < 2:
        raise ValueError("each study training partition requires at least two HC and two MCI rows")
    case_mean = np.nanmean(cases, axis=0)
    control_mean = np.nanmean(controls, axis=0)
    case_var = np.nanvar(cases, axis=0, ddof=1)
    control_var = np.nanvar(controls, axis=0, ddof=1)
    pooled_sd = np.sqrt(
        ((len(cases) - 1) * case_var + (len(controls) - 1) * control_var)
        / max(len(cases) + len(controls) - 2, 1)
    )
    effect = np.divide(
        case_mean - control_mean,
        pooled_sd,
        out=np.zeros(expression.shape[1], dtype=float),
        where=np.isfinite(pooled_sd) & (pooled_sd > 0),
    )
    effect[~np.isfinite(effect)] = 0.0
    return pd.Series(effect, index=expression.columns, name="effect")


def select_concordant_genes(
    effects: pd.DataFrame, top_k: int = TOP_K, allowed_genes: set[str] | None = None
) -> pd.DataFrame:
    """Select non-zero, same-direction study effects with deterministic stability ranking."""
    if list(effects.columns) != ["GSE63060", "GSE63061"]:
        raise ValueError("effects must contain GSE63060 and GSE63061 columns in that order")
    result = effects.copy()
    result.index.name = "gene"
    if allowed_genes is not None:
        result = result.loc[result.index.intersection(sorted(allowed_genes))].copy()
    result["direction_concordant"] = (
        np.isfinite(result["GSE63060"])
        & np.isfinite(result["GSE63061"])
        & (result["GSE63060"] * result["GSE63061"] > 0)
    )
    result = result.loc[result["direction_concordant"]].copy()
    result["direction"] = np.where(result["GSE63060"] > 0, "positive", "negative")
    result["stability_score"] = result[["GSE63060", "GSE63061"]].abs().min(axis=1)
    result["mean_abs_effect"] = result[["GSE63060", "GSE63061"]].abs().mean(axis=1)
    result = result.reset_index().sort_values(
        ["stability_score", "mean_abs_effect", "gene"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    result["selection_rank"] = np.arange(1, len(result) + 1)
    return result.head(top_k).reset_index(drop=True)


def stable_estimator(seed: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
        ]
    )


def pooled_estimator(feature_count: int, seed: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("select", SelectKBest(score_func=f_classif, k=min(TOP_K, feature_count))),
            ("model", LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=seed)),
        ]
    )


def development_endpoint() -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray, list[str]]:
    x60, meta60 = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    x61, meta61 = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    keep60 = meta60["label"].isin(LABELS)
    keep61 = meta61["label"].isin(LABELS)
    common = sorted(set(x60.columns).intersection(x61.columns))
    x60 = x60.loc[keep60, common]
    x61 = x61.loc[keep61, common]
    y60 = meta60.loc[x60.index, "label"].map(LABELS).to_numpy(dtype=int)
    y61 = meta61.loc[x61.index, "label"].map(LABELS).to_numpy(dtype=int)
    if (len(x60), int(y60.sum()), len(x61), int(y61.sum())) != (106, 39, 102, 30):
        raise ValueError("unexpected AddNeuroMed MCI/HC endpoint counts")
    return x60, y60, x61, y61, common


def selected_effects(
    x60: pd.DataFrame,
    y60: np.ndarray,
    x61: pd.DataFrame,
    y61: np.ndarray,
    allowed_genes: set[str],
) -> pd.DataFrame:
    return select_concordant_genes(
        pd.concat(
            [
                study_effects(x60, y60).rename("GSE63060"),
                study_effects(x61, y61).rename("GSE63061"),
            ],
            axis=1,
        ),
        allowed_genes=allowed_genes,
    )


def synchronized_outer_predictions(
    x60: pd.DataFrame,
    y60: np.ndarray,
    x61: pd.DataFrame,
    y61: np.ndarray,
    allowed_genes: set[str],
    seed: int,
    record_selection: bool,
    include_pooled: bool = True,
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    split60 = list(StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=seed).split(x60, y60))
    split61 = list(StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=seed).split(x61, y61))
    predictions = {
        "stable_60": np.full(len(y60), np.nan),
        "stable_61": np.full(len(y61), np.nan),
    }
    if include_pooled:
        predictions["pooled_60"] = np.full(len(y60), np.nan)
        predictions["pooled_61"] = np.full(len(y61), np.nan)
    events: list[dict[str, Any]] = []
    for fold, ((train60, test60), (train61, test61)) in enumerate(zip(split60, split61), start=1):
        selected = selected_effects(
            x60.iloc[train60], y60[train60], x61.iloc[train61], y61[train61], allowed_genes
        )
        genes = selected["gene"].tolist()
        if not genes:
            raise ValueError(f"fold {fold} has no direction-concordant genes")
        train_x = pd.concat([x60.iloc[train60][genes], x61.iloc[train61][genes]], axis=0)
        train_y = np.concatenate([y60[train60], y61[train61]])
        stable = stable_estimator(seed + fold)
        stable.fit(train_x, train_y)
        predictions["stable_60"][test60] = stable.predict_proba(x60.iloc[test60][genes])[:, 1]
        predictions["stable_61"][test61] = stable.predict_proba(x61.iloc[test61][genes])[:, 1]

        if include_pooled:
            model_genes = sorted(allowed_genes)
            pooled_x = pd.concat([x60.iloc[train60][model_genes], x61.iloc[train61][model_genes]], axis=0)
            pooled = pooled_estimator(pooled_x.shape[1], seed + 1000 + fold)
            pooled.fit(pooled_x, train_y)
            predictions["pooled_60"][test60] = pooled.predict_proba(x60.iloc[test60][model_genes])[:, 1]
            predictions["pooled_61"][test61] = pooled.predict_proba(x61.iloc[test61][model_genes])[:, 1]
        if record_selection:
            for record in selected.to_dict(orient="records"):
                events.append(
                    {
                        "fold": fold,
                        "training_gse63060_n": int(len(train60)),
                        "training_gse63061_n": int(len(train61)),
                        **record,
                    }
                )
    if any(np.isnan(values).any() for values in predictions.values()):
        raise ValueError("synchronized outer folds did not assign every held-out prediction")
    return predictions, events


def load_gse140829_unlabelled(allowed_genes: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Map external expression without reading any disease-bearing metadata fields."""
    probe_symbols = load_gpl_symbol_map(GPL15988_FULL)
    sample_ids: list[str] = []
    gene_sum: dict[str, np.ndarray] = {}
    gene_count: dict[str, np.ndarray] = {}
    in_table = False
    with gzip.open(GSE140829_MATRIX, "rt", encoding="utf-8", errors="replace", newline="") as handle:
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
            values = np.asarray([_numeric(value) for value in fields[1 : len(sample_ids) + 1]], dtype=np.float32)
            if len(values) != len(sample_ids):
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
    if len(expression) != 587:
        raise ValueError(f"expected 587 unlabelled GSE140829 rows, found {len(expression)}")
    return expression, {
        "series_sample_count": len(sample_ids),
        "mapped_requested_gene_count": len(gene_sum),
        "disease_metadata_read_before_prediction": False,
    }


def require_exact_features(expression: pd.DataFrame, genes: list[str], cohort: str) -> pd.DataFrame:
    missing = sorted(set(genes).difference(expression.columns))
    if missing:
        raise ValueError(f"{cohort} blocks complete locked scoring; missing {len(missing)} selected exact symbols")
    return expression.loc[:, genes]


def add_metric(
    y: np.ndarray,
    probability: np.ndarray,
    status: str,
    cohort: str,
    model: str,
    seed: int,
    **extra: Any,
) -> dict[str, Any]:
    row = metric_row(y, probability, status, "stage_mci_vs_hc", cohort, model)
    row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y, probability, seed)
    row.update(extra)
    return row


def main() -> int:
    x60, y60, x61, y61, common = development_endpoint()
    # This mapping-only candidate universe is fixed before any target diagnosis
    # metadata is read. It prevents a partial external score masquerading as a
    # complete stable-signal transport model.
    ext140_all, ext140_unlabelled_audit = load_gse140829_unlabelled(set(common))
    ext249_all, ext249_ids, ext249_unlabelled_audit = load_target_expression()
    transportable_genes = set(common).intersection(ext140_all.columns).intersection(ext249_all.columns)
    outer_predictions, selection_events = synchronized_outer_predictions(
        x60, y60, x61, y61, transportable_genes, SEED, True
    )
    metrics: list[dict[str, Any]] = []
    output_predictions: list[dict[str, Any]] = []
    for study, y, stable_key, pooled_key in [
        ("GSE63060", y60, "stable_60", "pooled_60"),
        ("GSE63061", y61, "stable_61", "pooled_61"),
    ]:
        for model, key in [("cross_study_stable_ridge", stable_key), ("pooled_ridge_nested_selector", pooled_key)]:
            probability = outer_predictions[key]
            metrics.append(
                add_metric(
                    y,
                    probability,
                    "DEVELOPMENT_FAMILY_SYNCHRONIZED_OUTER_FOLD",
                    study,
                    model,
                    SEED + len(metrics),
                    evidence_tier="development-family cross-study consistency; not independent validation",
                    target_labels_used_for_selection=False,
                    target_labels_used_for_fit=False,
                )
            )
            output_predictions.extend(
                {
                    "cohort": study,
                    "assessment_index": int(index),
                    "y": int(label),
                    "probability": float(probability[index]),
                    "model": model,
                    "status": "DEVELOPMENT_FAMILY_SYNCHRONIZED_OUTER_FOLD",
                }
                for index, label in enumerate(y)
            )

    final_selected = selected_effects(x60, y60, x61, y61, transportable_genes)
    final_genes = final_selected["gene"].tolist()
    full_x = pd.concat([x60[final_genes], x61[final_genes]], axis=0)
    full_y = np.concatenate([y60, y61])
    stable_final = stable_estimator(SEED + 999)
    stable_final.fit(full_x, full_y)
    model_genes = sorted(transportable_genes)
    pooled_full_x = pd.concat([x60[model_genes], x61[model_genes]], axis=0)
    pooled_final = pooled_estimator(pooled_full_x.shape[1], SEED + 1999)
    pooled_final.fit(pooled_full_x, full_y)

    # Both target matrices were loaded without diagnosis metadata before selection.
    stable_140_all = stable_final.predict_proba(require_exact_features(ext140_all, final_genes, "GSE140829"))[:, 1]
    pooled_140_all = pooled_final.predict_proba(require_exact_features(ext140_all, model_genes, "GSE140829"))[:, 1]
    stable_249_all = stable_final.predict_proba(require_exact_features(ext249_all, final_genes, "GSE249477"))[:, 1]
    pooled_249_all = pooled_final.predict_proba(require_exact_features(ext249_all, model_genes, "GSE249477"))[:, 1]

    # Assessment metadata are now read only after all four target prediction vectors are fixed.
    _, ext140_meta, ext140_assessment_audit = load_external(
        GSE140829_MATRIX, GPL15988_FULL, set(common), include_unlabelled_rows=True
    )
    if ext140_meta.index.tolist() != ext140_all.index.tolist():
        raise ValueError("GSE140829 assessment metadata do not preserve preassessment row order")
    labels140 = ext140_meta["label"]
    labels249 = parse_gse249477_assessment_labels(ext249_ids)
    for cohort, labels, probabilities, all_rows, expected_n in [
        ("GSE140829", labels140, {"cross_study_stable_ridge": stable_140_all, "pooled_ridge_nested_selector": pooled_140_all}, 587, 383),
        ("GSE249477", labels249, {"cross_study_stable_ridge": stable_249_all, "pooled_ridge_nested_selector": pooled_249_all}, 62, 41),
    ]:
        keep = labels.isin(LABELS)
        y = labels.loc[keep].map(LABELS).to_numpy(dtype=int)
        if len(y) != expected_n:
            raise ValueError(f"{cohort} assessment endpoint count changed: {len(y)}")
        for model, probability_all in probabilities.items():
            probability = probability_all[np.asarray(keep, dtype=bool)]
            metrics.append(
                add_metric(
                    y,
                    probability,
                    "EXTERNAL_TRANSPORT_AUDIT" if cohort == "GSE140829" else "SECONDARY_INDEPENDENT_CROSS_PLATFORM_TRANSPORT_AUDIT",
                    cohort,
                    model,
                    SEED + len(metrics),
                    all_rows_predicted_before_endpoint_assessment=all_rows,
                    target_labels_used_for_selection=False,
                    target_labels_used_for_fit=False,
                    target_labels_used_for_preprocessing=False,
                    target_labels_used_for_threshold_or_prediction_calibration=False,
                    target_labels_used_for_descriptive_assessment=True,
                    target_labels_used_only_for_assessment=True,
                    candidate_gene_universe="mapping_only_exact_symbol_intersection",
                    evidence_tier="locked external transport audit; not external clinical validation",
                )
            )
            output_predictions.extend(
                {
                    "cohort": cohort,
                    "assessment_index": int(index),
                    "y": int(label),
                    "probability": float(probability[index]),
                    "model": model,
                    "status": "EXTERNAL_TRANSPORT_AUDIT" if cohort == "GSE140829" else "SECONDARY_INDEPENDENT_CROSS_PLATFORM_TRANSPORT_AUDIT",
                }
                for index, label in enumerate(y)
            )

    observed = float(
        add_metric(
            np.concatenate([y60, y61]),
            np.concatenate([outer_predictions["stable_60"], outer_predictions["stable_61"]]),
            "DEVELOPMENT_FAMILY_SYNCHRONIZED_OUTER_FOLD",
            "GSE63060_GSE63061",
            "cross_study_stable_ridge",
            SEED + 700,
        )["auc"]
    )
    null_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(SEED + 5000)
    for repeat in range(1, PERMUTATION_REPETITIONS + 1):
        permuted60 = rng.permutation(y60)
        permuted61 = rng.permutation(y61)
        permuted, _ = synchronized_outer_predictions(
            x60, permuted60, x61, permuted61, transportable_genes, SEED, False, include_pooled=False
        )
        null_auc = add_metric(
            np.concatenate([permuted60, permuted61]),
            np.concatenate([permuted["stable_60"], permuted["stable_61"]]),
            "DEVELOPMENT_LABEL_PERMUTATION_NULL",
            "GSE63060_GSE63061",
            "cross_study_stable_ridge",
            SEED + 8000 + repeat,
        )["auc"]
        null_rows.append({"repeat": repeat, "development_oof_auc": null_auc})
    null = pd.DataFrame(null_rows)
    permutation_p = float((1 + int((null["development_oof_auc"] >= observed).sum())) / (1 + len(null)))

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metrics).to_csv(OUT / "cross_study_stable_metrics.csv", index=False)
    pd.DataFrame(output_predictions).to_csv(OUT / "cross_study_stable_predictions_anonymized.csv", index=False)
    final_selected.assign(selection_scope="all eligible AddNeuroMed MCI/HC rows").to_csv(
        OUT / "cross_study_stable_selected_genes.csv", index=False
    )
    pd.DataFrame(selection_events).to_csv(OUT / "cross_study_stable_selection_events.csv", index=False)
    null.to_csv(OUT / "cross_study_stable_permutation_null.csv", index=False)
    audit = {
        "status": "CROSS_STUDY_STABLE_SIGNAL_AUDIT",
        "design": "synchronized development-family outer folds; direction-concordant effect selection in paired training partitions only",
        "endpoint": "stage_mci_vs_hc",
        "development": {"GSE63060": {"rows": len(y60), "HC": int((y60 == 0).sum()), "MCI": int(y60.sum())}, "GSE63061": {"rows": len(y61), "HC": int((y61 == 0).sum()), "MCI": int(y61.sum())}},
        "final_selected_gene_count": len(final_genes),
        "mapping_only_transportable_gene_universe_count": len(transportable_genes),
        "final_selection_direction_concordant": bool(final_selected["direction_concordant"].all()),
        "external_assessment": {"GSE140829_rows": 383, "GSE249477_rows": 41},
        "gse140829_all_rows_predicted_before_endpoint_assessment": 587,
        "gse249477_all_rows_predicted_before_endpoint_assessment": 62,
        "gse140829_target_labels_used_for_selection": False,
        "gse140829_target_labels_used_for_fit": False,
        "gse140829_target_labels_used_for_preprocessing": False,
        "gse140829_target_labels_used_for_descriptive_assessment": True,
        "gse140829_target_labels_used_only_for_assessment": True,
        "gse249477_target_labels_used_for_selection": False,
        "gse249477_target_labels_used_for_fit": False,
        "gse249477_target_labels_used_for_preprocessing": False,
        "gse249477_target_labels_used_for_descriptive_assessment": True,
        "gse249477_target_labels_used_only_for_assessment": True,
        "prediction_calibration_or_adaptation_fit_on_target": False,
        "permutation_null_repetitions": PERMUTATION_REPETITIONS,
        "development_oof_auc": observed,
        "development_label_permutation_p": permutation_p,
        "gse140829_unlabelled_expression_audit": ext140_unlabelled_audit,
        "gse140829_assessment_metadata_audit": ext140_assessment_audit,
        "gse249477_unlabelled_expression_audit": ext249_unlabelled_audit,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [
                GSE60_MATRIX,
                GSE60_LABELS,
                GSE61_MATRIX,
                GSE61_LABELS,
                GSE140829_MATRIX,
                GPL15988_FULL,
                GSE249477_RAW,
                GSE249477_SERIES,
                Path(__file__),
            ]
        },
        "allowed_claim": "direction-concordant development-family signals were stress-tested under locked external transport assessment",
        "forbidden_claim": "the audit externally validates a clinically deployable MCI biomarker or solves domain shift",
    }
    (OUT / "cross_study_stable_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"metrics": len(metrics), "selected_genes": len(final_genes), "permutation_p": permutation_p}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
