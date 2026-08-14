from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import f_classif
from sklearn.model_selection import StratifiedKFold

from run_gse249477_rank_transport import (
    GSE249477_RAW,
    GSE249477_SERIES,
    load_target_expression,
    parse_assessment_labels,
)
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
SEED = 20260813
CANDIDATE_GENE_COUNT = 100
PAIR_COUNT = 10


def select_tsp_pairs(
    expression: pd.DataFrame,
    y: np.ndarray,
    candidate_gene_count: int = CANDIDATE_GENE_COUNT,
    pair_count: int = PAIR_COUNT,
) -> list[tuple[str, str, int]]:
    """Select deterministic, disjoint MCI-oriented relative-expression pairs."""
    if len(np.unique(y)) != 2:
        raise ValueError("TSP selection requires two development classes")
    values = expression.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("TSP selection requires finite development expression")
    scores, _ = f_classif(values, y)
    candidates = [
        gene
        for gene, score in sorted(
            zip(expression.columns.astype(str), scores), key=lambda item: (-np.nan_to_num(item[1], nan=-np.inf), item[0])
        )
        if np.isfinite(score)
    ][: min(candidate_gene_count, expression.shape[1])]
    if len(candidates) < 2:
        raise ValueError("TSP selection needs at least two finite-score genes")

    mci = expression.iloc[np.flatnonzero(y == 1)]
    hc = expression.iloc[np.flatnonzero(y == 0)]
    pair_rows: list[tuple[float, str, str, int]] = []
    for left, right in combinations(candidates, 2):
        gap = float((mci[left] > mci[right]).mean() - (hc[left] > hc[right]).mean())
        if gap > 0:
            pair_rows.append((abs(gap), left, right, 1))
        elif gap < 0:
            pair_rows.append((abs(gap), left, right, -1))
    pair_rows.sort(key=lambda row: (-row[0], row[1], row[2], row[3]))
    selected: list[tuple[str, str, int]] = []
    used: set[str] = set()
    for _gap, left, right, orientation in pair_rows:
        if left in used or right in used:
            continue
        selected.append((left, right, orientation))
        used.update({left, right})
        if len(selected) == pair_count:
            break
    if len(selected) != pair_count:
        raise ValueError(f"expected {pair_count} disjoint TSP pairs, selected {len(selected)}")
    return selected


def predict_pair_vote(expression: pd.DataFrame, pairs: list[tuple[str, str, int]]) -> np.ndarray:
    """Score frozen relative-expression pairs without labels or target fitting."""
    if not pairs:
        raise ValueError("TSP prediction requires at least one pair")
    votes = []
    for left, right, orientation in pairs:
        if left not in expression.columns or right not in expression.columns:
            raise ValueError(f"TSP target expression lacks frozen pair: {left}, {right}")
        ordering = expression[left].to_numpy(dtype=float) > expression[right].to_numpy(dtype=float)
        votes.append(ordering if orientation == 1 else ~ordering)
    return np.mean(np.vstack(votes), axis=0, dtype=float)


def load_development() -> tuple[pd.DataFrame, np.ndarray]:
    dev60_x, dev60_meta = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    dev61_x, dev61_meta = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    common = sorted(set(dev60_x.columns).intersection(dev61_x.columns))
    expression = pd.concat([dev60_x[common], dev61_x[common]], axis=0)
    metadata = pd.concat([dev60_meta, dev61_meta], axis=0).loc[expression.index]
    keep = metadata["label"].isin({"HC", "MCI"})
    expression = expression.loc[keep].astype(np.float32)
    if expression.isna().any().any():
        expression = expression.fillna(expression.median(axis=0)).fillna(0.0)
    return expression, metadata.loc[keep, "label"].map({"HC": 0, "MCI": 1}).to_numpy(dtype=int)


def assess(
    y: np.ndarray,
    probability: np.ndarray,
    status: str,
    cohort: str,
    seed: int,
) -> dict[str, object]:
    row = metric_row(y, probability, status, "stage_mci_vs_hc", cohort, "relative_expression_tsp")
    row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y, probability, seed)
    row["candidate_gene_count"] = CANDIDATE_GENE_COUNT
    row["pair_count"] = PAIR_COUNT
    if cohort == "GSE249477":
        row["calibration_intercept"] = np.nan
        row["calibration_slope"] = np.nan
        row["calibration_method"] = "NOT_FIT_TARGET_LABELS"
        row["calibration_fit_scope"] = "not_fit_target_labels"
    return row


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path == Path(__file__).resolve() else sha256_file(path)


def main(output_dir: Path = OUT) -> int:
    dev_expression, y_dev = load_development()
    target_expression, _target_ids, target_input_audit = load_target_expression()
    external_expression, external_meta, external_input_audit = load_external(
        GSE140829_MATRIX, GPL15988_FULL, set(dev_expression.columns)
    )
    common = sorted(set(dev_expression.columns).intersection(external_expression.columns).intersection(target_expression.columns))
    if len(common) < 2 * PAIR_COUNT:
        raise ValueError(f"need at least {2 * PAIR_COUNT} shared genes, found {len(common)}")
    dev_expression = dev_expression[common]
    external_expression = external_expression[common].fillna(0.0)
    target_expression = target_expression[common].fillna(0.0)

    folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = np.full(len(y_dev), np.nan, dtype=float)
    fold_pairs: list[dict[str, object]] = []
    for fold, (train_idx, test_idx) in enumerate(folds.split(dev_expression, y_dev), start=1):
        pairs = select_tsp_pairs(dev_expression.iloc[train_idx], y_dev[train_idx])
        oof[test_idx] = predict_pair_vote(dev_expression.iloc[test_idx], pairs)
        fold_pairs.extend(
            {"selection_scope": "development_fold", "fold": fold, "pair_rank": rank, "gene_a": left, "gene_b": right, "mci_orientation": orientation}
            for rank, (left, right, orientation) in enumerate(pairs, start=1)
        )
    final_pairs = select_tsp_pairs(dev_expression, y_dev)

    # GSE140829 expression is scored before restricting to the endpoint labels.
    external_probability_all = predict_pair_vote(external_expression, final_pairs)
    external_mask = external_meta["label"].isin({"HC", "MCI"})
    y_external = external_meta.loc[external_mask, "label"].map({"HC": 0, "MCI": 1}).to_numpy(dtype=int)
    external_probability = external_probability_all[external_mask.to_numpy()]

    # GSE249477 expression is scored across every target row before diagnosis parsing.
    target_probability_all = predict_pair_vote(target_expression, final_pairs)
    target_labels = parse_assessment_labels(target_expression.index.tolist())
    target_mask = target_labels.isin({"HC", "MCI"})
    y_target = target_labels.loc[target_mask].map({"HC": 0, "MCI": 1}).to_numpy(dtype=int)
    target_probability = target_probability_all[target_mask.to_numpy()]

    internal = assess(y_dev, oof, "NESTED_INTERNAL_SENSITIVITY", "GSE63060_GSE63061", SEED + 100)
    external = assess(y_external, external_probability, "PRIMARY_EXTERNAL_TRANSPORT_AUDIT", "GSE140829", SEED + 200)
    secondary = assess(
        y_target,
        target_probability,
        "SECONDARY_INDEPENDENT_CROSS_PLATFORM_TRANSPORT_AUDIT",
        "GSE249477",
        SEED + 300,
    )
    for row in [external, secondary]:
        row["target_labels_used_for_pair_selection"] = False
        row["target_labels_used_for_prediction"] = False
        row["target_labels_used_for_threshold_or_calibration"] = False
        row["target_labels_used_only_for_assessment"] = True

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([internal, external, secondary]).to_csv(output_dir / "tsp_transport_metrics.csv", index=False)
    pair_rows = fold_pairs + [
        {
            "selection_scope": "final_development_fit",
            "fold": "all",
            "pair_rank": rank,
            "gene_a": left,
            "gene_b": right,
            "mci_orientation": orientation,
        }
        for rank, (left, right, orientation) in enumerate(final_pairs, start=1)
    ]
    pd.DataFrame(pair_rows).to_csv(output_dir / "tsp_transport_selected_pairs.csv", index=False)
    pd.concat(
        [
            pd.DataFrame({"cohort": "GSE140829", "assessment_index": np.arange(len(y_external)), "y": y_external, "probability": external_probability}),
            pd.DataFrame({"cohort": "GSE249477", "assessment_index": np.arange(len(y_target)), "y": y_target, "probability": target_probability}),
        ],
        ignore_index=True,
    ).to_csv(output_dir / "tsp_transport_predictions_anonymized.csv", index=False)
    audit = {
        "status": "RELATIVE_EXPRESSION_TRANSPORT_AUDIT",
        "method": "development-only Top-Scoring Pairs relative-expression vote",
        "fixed_candidate_gene_count": CANDIDATE_GENE_COUNT,
        "fixed_pair_count": PAIR_COUNT,
        "development": {"cohort": "GSE63060_GSE63061", "rows": int(len(y_dev)), "mci": int(y_dev.sum()), "hc": int((y_dev == 0).sum())},
        "gse140829": {"assessment_rows": int(len(y_external)), "label_counts": {"HC": int((y_external == 0).sum()), "MCI": int(y_external.sum())}, "all_rows_predicted_before_endpoint_filter": int(len(external_probability_all))},
        "gse249477": {"assessment_rows": int(len(y_target)), "label_counts": {"HC": int((y_target == 0).sum()), "MCI": int(y_target.sum())}, "all_rows_predicted_before_label_parsing": int(len(target_probability_all)), "excluded_ad_rows": int((target_labels == "AD").sum())},
        "gse140829_target_labels_used_for_pair_selection": False,
        "gse140829_target_labels_used_for_prediction": False,
        "gse140829_target_labels_used_for_threshold_or_calibration": False,
        "gse140829_target_labels_used_only_for_assessment": True,
        "gse249477_target_labels_used_for_pair_selection": False,
        "gse249477_target_labels_used_for_prediction": False,
        "gse249477_target_labels_used_for_threshold_or_calibration": False,
        "gse249477_target_labels_used_only_for_assessment": True,
        "gse249477_identifier_audit": target_input_audit,
        "gse140829_input_audit": external_input_audit,
        "shared_exact_case_canonicalized_gene_count": int(len(common)),
        "final_pairs_are_disjoint": len({gene for pair in final_pairs for gene in pair[:2]}) == 2 * PAIR_COUNT,
        "allowed_claims": [
            "A development-only relative-expression TSP method was tested as a transport audit.",
            "Both target domains were scored before their endpoint labels were used for aggregate assessment.",
        ],
        "forbidden_claims": [
            "TSP externally validates an MCI clinical biomarker.",
            "Relative expression ordering solves domain shift.",
            "TSP probabilities are target calibrated or clinically deployable.",
        ],
        "source_hashes": {
            str(path): source_hash(path)
            for path in [
                GSE60_MATRIX,
                GSE60_LABELS,
                GSE61_MATRIX,
                GSE61_LABELS,
                GSE140829_MATRIX,
                GPL15988_FULL,
                GSE249477_RAW,
                GSE249477_SERIES,
                Path(__file__).resolve(),
            ]
        },
    }
    (output_dir / "tsp_transport_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(pd.DataFrame([internal, external, secondary]).to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    raise SystemExit(main(args.output_dir))
