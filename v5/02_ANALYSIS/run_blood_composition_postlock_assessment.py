"""Post-lock exploratory MCI-versus-HC assessment of frozen MCPcounter scores."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests

from run_blood_composition_audit import (
    GSE140829_MATRIX,
    GSE63060_MATRIX,
    GSE63061_MATRIX,
    GSE249477_RAW,
    GPL15988_FULL,
    ROOT,
    sha256_file,
)
from run_gse249477_rank_transport import parse_assessment_labels
from run_gse249477_rank_transport import GSE249477_SERIES
from run_pairwise_external import (
    GSE60_LABELS,
    GSE61_LABELS,
    load_dev,
    load_external,
)


OUT = ROOT / "03_RESULTS"
PREASSESSMENT_SCORES = OUT / "blood_composition_preassessment_scores_anonymized.csv"
PREASSESSMENT_MANIFEST = OUT / "blood_composition_preassessment_manifest.json"
PROTOCOL = ROOT / "00_PROTOCOL" / "blood_composition_postlock_assessment_v1.md"
POSTLOCK_SOURCE_REGISTRY = ROOT / "00_PROTOCOL" / "blood_composition_postlock_source_registry_v1.json"
SUMMARY_NAME = "blood_composition_postlock_summary.csv"
CORRELATIONS_NAME = "blood_composition_postlock_prediction_correlations.csv"
AUDIT_NAME = "blood_composition_postlock_audit.json"
LABELS = {"HC": 0, "MCI": 1}
COHORTS = ("GSE63060", "GSE63061", "GSE140829", "GSE249477")
LOCKED_PREDICTIONS = OUT / "cross_study_stable_predictions_anonymized.csv"
LOCKED_PREDICTION_AUDIT = OUT / "cross_study_stable_audit.json"
LABEL_SOURCE_PATHS = {
    "GSE63060": GSE60_LABELS,
    "GSE63061": GSE61_LABELS,
    "GSE140829": GSE140829_MATRIX,
    "GSE249477": GSE249477_SERIES,
}
EXPRESSION_SOURCE_PATHS = {
    "GSE63060": GSE63060_MATRIX,
    "GSE63061": GSE63061_MATRIX,
    "GSE140829": GSE140829_MATRIX,
    "GSE140829_platform": GPL15988_FULL,
    "GSE249477": GSE249477_RAW,
}


def hedges_g(mci: np.ndarray, hc: np.ndarray) -> float:
    """Return the small-sample-corrected MCI-minus-HC standardized mean difference."""
    n_mci = len(mci)
    n_hc = len(hc)
    degrees_of_freedom = n_mci + n_hc - 2
    if n_mci < 2 or n_hc < 2 or degrees_of_freedom <= 0:
        return np.nan
    mci_variance = np.var(mci, ddof=1)
    hc_variance = np.var(hc, ddof=1)
    pooled_sd = np.sqrt(
        ((n_mci - 1) * mci_variance + (n_hc - 1) * hc_variance) / degrees_of_freedom
    )
    if not np.isfinite(pooled_sd) or pooled_sd == 0:
        return 0.0 if np.isclose(np.mean(mci), np.mean(hc)) else np.nan
    cohen_d = (np.mean(mci) - np.mean(hc)) / pooled_sd
    correction = 1.0 - 3.0 / (4.0 * degrees_of_freedom - 1.0)
    return float(correction * cohen_d)


def _direction(mean_difference: float) -> str:
    if np.isclose(mean_difference, 0.0):
        return "no_mean_difference"
    return "MCI_higher" if mean_difference > 0 else "HC_higher"


def _paths_match(expected: Path, registered: object) -> bool:
    """Compare registry paths without resolving a nonexistent synthetic test path."""
    return resolve_registry_path(str(registered)) == Path(expected).absolute()


def resolve_registry_path(path_text: str) -> Path:
    """Resolve registry-relative paths from the project root, not the caller's CWD."""
    path = Path(path_text)
    return path.absolute() if path.is_absolute() else (ROOT / path).absolute()


def _validate_registry_entry(
    entry: object,
    expected_path: Path,
    purpose: str,
) -> dict[str, Any]:
    if not isinstance(entry, dict) or not isinstance(entry.get("path"), str) or not isinstance(entry.get("sha256"), str):
        raise ValueError(f"{purpose} registry entry lacks path and SHA256")
    if not _paths_match(expected_path, entry["path"]):
        raise ValueError(f"{purpose} registry path does not match the declared source")
    if not expected_path.exists():
        raise FileNotFoundError(f"{purpose} source is unavailable: {expected_path}")
    if sha256_file(expected_path).lower() != entry["sha256"].lower():
        raise ValueError(f"{purpose} SHA256 does not match the frozen registry")
    return entry


def load_postlock_source_registry() -> dict[str, Any]:
    """Load the versioned source registry before any post-lock label access."""
    if not POSTLOCK_SOURCE_REGISTRY.exists():
        raise FileNotFoundError("post-lock source registry is required before diagnosis assessment")
    registry = json.loads(POSTLOCK_SOURCE_REGISTRY.read_text(encoding="utf-8"))
    if registry.get("status") != "POSTLOCK_COMPOSITION_SOURCE_REGISTRY_V1":
        raise ValueError("post-lock source registry has an unexpected status")
    return registry


def validate_label_sources(registry: dict[str, Any]) -> None:
    """Fail closed on any diagnosis-source substitution before parsing labels."""
    entries = registry.get("label_sources")
    if not isinstance(entries, dict) or set(entries) != set(COHORTS):
        raise ValueError("post-lock source registry lacks the complete label_sources contract")
    for cohort in COHORTS:
        _validate_registry_entry(entries[cohort], LABEL_SOURCE_PATHS[cohort], f"{cohort} label source")


def validate_expression_sources(preassessment_manifest: dict[str, Any]) -> None:
    """Confirm the all-row score-to-label row mappings still use the frozen expression inputs."""
    input_hashes = preassessment_manifest.get("input_hashes")
    if not isinstance(input_hashes, dict):
        raise ValueError("preassessment manifest lacks expression input hashes")
    for source_name, source_path in EXPRESSION_SOURCE_PATHS.items():
        expected_hash = input_hashes.get(str(source_path))
        if not isinstance(expected_hash, str):
            raise ValueError(f"preassessment manifest lacks {source_name} expression-source SHA256")
        if not source_path.exists():
            raise FileNotFoundError(f"{source_name} expression source is unavailable: {source_path}")
        if sha256_file(source_path).lower() != expected_hash.lower():
            raise ValueError(f"{source_name} expression source SHA256 does not match the frozen preassessment manifest")


def load_locked_predictions(registry: dict[str, Any]) -> pd.DataFrame:
    """Return immutable prediction rows only after prediction and audit provenance checks."""
    artifacts = registry.get("prediction_artifacts")
    if not isinstance(artifacts, dict):
        raise ValueError("post-lock source registry lacks prediction_artifacts")
    _validate_registry_entry(artifacts.get("predictions"), LOCKED_PREDICTIONS, "locked prediction")
    audit_entry = _validate_registry_entry(
        artifacts.get("audit"), LOCKED_PREDICTION_AUDIT, "locked prediction audit"
    )
    prediction_audit = json.loads(LOCKED_PREDICTION_AUDIT.read_text(encoding="utf-8"))
    if prediction_audit.get("status") != audit_entry.get("required_status"):
        raise ValueError("locked prediction audit status does not match the frozen registry")
    return pd.read_csv(LOCKED_PREDICTIONS)


def summarize_mci_hc_scores(
    cohort: str,
    scores: pd.DataFrame,
    labels: pd.Series,
) -> pd.DataFrame:
    """Summarize every frozen score with a within-cohort MCPcounter FDR family."""
    if not scores.index.equals(labels.index):
        raise ValueError(f"{cohort} score and diagnosis indices do not match exactly")
    endpoint = labels.isin(LABELS)
    if not endpoint.any():
        raise ValueError(f"{cohort} has no MCI/HC rows for post-lock assessment")
    endpoint_scores = scores.loc[endpoint]
    endpoint_labels = labels.loc[endpoint]
    mci_mask = endpoint_labels.eq("MCI").to_numpy()
    hc_mask = endpoint_labels.eq("HC").to_numpy()
    rows: list[dict[str, object]] = []
    for score_name in endpoint_scores.columns:
        values = endpoint_scores[score_name].to_numpy(dtype=float)
        mci = values[mci_mask]
        hc = values[hc_mask]
        if not np.isfinite(values).all():
            raise ValueError(f"{cohort} {score_name} contains a non-finite locked score")
        test = ttest_ind(mci, hc, equal_var=False, alternative="two-sided")
        mean_difference = float(np.mean(mci) - np.mean(hc))
        rows.append(
            {
                "cohort": cohort,
                "method": "MCPcounter",
                "cell_score": score_name,
                "analysis_unit": "public-data sample row",
                "n_hc": int(len(hc)),
                "n_mci": int(len(mci)),
                "mean_hc": float(np.mean(hc)),
                "mean_mci": float(np.mean(mci)),
                "mean_difference_mci_minus_hc": mean_difference,
                "hedges_g": hedges_g(mci, hc),
                "welch_t": float(test.statistic),
                "p_value_two_sided": float(test.pvalue),
                "direction": _direction(mean_difference),
            }
        )
    result = pd.DataFrame(rows)
    result["p_value_bh"] = multipletests(result["p_value_two_sided"], method="fdr_bh")[1]
    result["fdr_family"] = f"{cohort}|MCPcounter|{len(result)}_scores"
    result["exploratory_only"] = True
    return result.sort_values("cell_score", kind="mergesort").reset_index(drop=True)


def align_locked_predictions(
    score_endpoint: pd.DataFrame,
    predictions: pd.DataFrame,
    cohort: str,
    model: str,
) -> pd.DataFrame:
    """Join a pre-existing prediction only after a complete label/order agreement check."""
    required_scores = {"assessment_index", "y", "cell_score"}
    required_predictions = {"cohort", "assessment_index", "y", "probability", "model"}
    if not required_scores.issubset(score_endpoint.columns):
        raise ValueError("score endpoint table lacks the required locked-mapping columns")
    if not required_predictions.issubset(predictions.columns):
        raise ValueError("prediction table lacks the required locked-mapping columns")
    selected = predictions.loc[(predictions["cohort"] == cohort) & (predictions["model"] == model)].copy()
    if selected["assessment_index"].duplicated().any():
        raise ValueError(f"{cohort} {model} has duplicate locked prediction indices")
    expected = score_endpoint.loc[:, ["assessment_index", "y", "cell_score"]].copy()
    if expected["assessment_index"].duplicated().any():
        raise ValueError(f"{cohort} endpoint mapping has duplicate assessment indices")
    joined = expected.merge(
        selected.loc[:, ["assessment_index", "y", "probability"]],
        on="assessment_index",
        how="outer",
        suffixes=("_scores", "_prediction"),
        validate="one_to_one",
        indicator=True,
    )
    if not (joined["_merge"] == "both").all():
        raise ValueError(f"{cohort} {model} locked predictions do not cover the score endpoint exactly")
    if not (joined["y_scores"] == joined["y_prediction"]).all():
        raise ValueError(f"{cohort} {model} locked prediction labels do not match the score endpoint")
    return joined.loc[:, ["assessment_index", "y_scores", "cell_score", "probability"]].rename(
        columns={"y_scores": "y"}
    )


def _load_preassessment() -> tuple[pd.DataFrame, dict[str, Any]]:
    if not PREASSESSMENT_SCORES.exists() or not PREASSESSMENT_MANIFEST.exists():
        raise FileNotFoundError("locked MCPcounter preassessment artifacts are required before post-lock assessment")
    manifest = json.loads(PREASSESSMENT_MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "MCP_COUNTER_SOURCE_BRIDGE_PREASSESSMENT":
        raise ValueError("post-lock assessment requires a locked MCPcounter preassessment manifest")
    if manifest.get("diagnosis_labels_accessed") is not False:
        raise ValueError("preassessment manifest does not preserve diagnosis-label isolation")
    if manifest.get("method") != "MCPcounter" or manifest.get("xcell", {}).get("status") != "BLOCKED":
        raise ValueError("post-lock assessment requires the prespecified single-method MCPcounter state")
    observed_hash = sha256_file(PREASSESSMENT_SCORES)
    if observed_hash != manifest.get("anonymized_output", {}).get("sha256"):
        raise ValueError("locked preassessment score SHA256 does not match its manifest")
    scores = pd.read_csv(PREASSESSMENT_SCORES)
    validate_preassessment_contract(scores, manifest)
    return scores, manifest


def validate_preassessment_contract(scores: pd.DataFrame, manifest: dict[str, Any]) -> None:
    """Verify all score-table structural locks before any diagnosis source can be opened."""
    output = manifest.get("anonymized_output")
    cohort_preassessment = manifest.get("cohort_preassessment")
    if not isinstance(output, dict) or not isinstance(cohort_preassessment, dict):
        raise ValueError("preassessment manifest lacks anonymized-output and cohort_preassessment contracts")
    if output.get("path") != PREASSESSMENT_SCORES.name:
        raise ValueError("preassessment manifest anonymized-output path does not match the score table")
    if scores.columns.tolist() != output.get("columns"):
        raise ValueError("preassessment manifest columns do not match the locked score table")
    if len(scores) != output.get("row_count"):
        raise ValueError("preassessment manifest row count does not match the locked score table")
    score_columns = scores.columns[2:].tolist()
    if scores.columns[:2].tolist() != ["cohort", "assessment_index"] or len(score_columns) != 10:
        raise ValueError("locked MCPcounter score table has an unexpected column contract")
    if set(cohort_preassessment) != set(COHORTS):
        raise ValueError("preassessment manifest cohort_preassessment contract is incomplete")
    observed_cohorts = set(scores["cohort"])
    if observed_cohorts != set(COHORTS):
        raise ValueError("preassessment manifest cohort_preassessment does not match the score-table cohorts")
    for cohort in COHORTS:
        entry = cohort_preassessment[cohort]
        if not isinstance(entry, dict):
            raise ValueError(f"preassessment manifest cohort_preassessment entry is invalid for {cohort}")
        cohort_scores = scores.loc[scores["cohort"] == cohort]
        row_count = len(cohort_scores)
        if entry.get("expected_rows") != row_count or entry.get("scored_rows") != row_count:
            raise ValueError(f"preassessment manifest cohort_preassessment row coverage does not match {cohort}")
        if entry.get("score_column_count") != len(score_columns):
            raise ValueError(f"preassessment manifest cohort_preassessment score-column count does not match {cohort}")
        if cohort_scores["assessment_index"].tolist() != list(range(row_count)):
            raise ValueError(f"locked score table assessment indices are not contiguous and ordered for {cohort}")
    if not np.isfinite(scores.loc[:, score_columns].to_numpy(dtype=float)).all():
        raise ValueError("locked MCPcounter score table contains non-finite values")


def _label_series_by_preassessment_order(cohort: str, expected_rows: int) -> pd.Series:
    if cohort == "GSE63060":
        from run_blood_composition_audit import load_development_expression_label_blind

        expression = load_development_expression_label_blind(GSE63060_MATRIX, expected_rows, cohort)
        _development_expression, metadata = load_dev(GSE63060_MATRIX, GSE60_LABELS, cohort)
        return metadata["label"].reindex(expression.index)
    if cohort == "GSE63061":
        from run_blood_composition_audit import load_development_expression_label_blind

        expression = load_development_expression_label_blind(GSE63061_MATRIX, expected_rows, cohort)
        _development_expression, metadata = load_dev(GSE63061_MATRIX, GSE61_LABELS, cohort)
        return metadata["label"].reindex(expression.index)
    if cohort == "GSE140829":
        from run_blood_composition_audit import load_gse140829_expression_label_blind

        expression = load_gse140829_expression_label_blind()
        _external_expression, metadata, _audit = load_external(
            GSE140829_MATRIX,
            GPL15988_FULL,
            set(expression.columns),
            include_unlabelled_rows=True,
        )
        return metadata["label"].reindex(expression.index)
    if cohort == "GSE249477":
        from run_blood_composition_audit import load_gse249477_expression_label_blind

        expression, _audit = load_gse249477_expression_label_blind()
        return parse_assessment_labels(expression.index.tolist())
    raise ValueError(f"unsupported post-lock cohort: {cohort}")


def _cohort_scores_with_labels(all_scores: pd.DataFrame, cohort: str) -> tuple[pd.DataFrame, pd.Series]:
    subset = all_scores.loc[all_scores["cohort"] == cohort].copy()
    if subset.empty:
        raise ValueError(f"locked score table lacks cohort {cohort}")
    expected_rows = len(subset)
    expected_index = list(range(expected_rows))
    if subset["assessment_index"].tolist() != expected_index:
        raise ValueError(f"{cohort} locked score table assessment indices are not contiguous and ordered")
    labels = _label_series_by_preassessment_order(cohort, expected_rows)
    if len(labels) != expected_rows:
        raise ValueError(f"{cohort} diagnosis sequence does not cover locked preassessment rows")
    labels = pd.Series(labels.to_numpy(), index=subset["assessment_index"].to_numpy(), name="label")
    score_frame = subset.drop(columns=["cohort", "assessment_index"])
    score_frame.index = labels.index
    return score_frame, labels


def _endpoint_mapping(labels: pd.Series) -> pd.DataFrame:
    endpoint = labels.isin(LABELS)
    return pd.DataFrame(
        {
            "assessment_index": labels.index[endpoint].to_numpy(dtype=int),
            "y": labels.loc[endpoint].map(LABELS).to_numpy(dtype=int),
        }
    )


def _describe_prediction_correlations(
    cohort: str,
    endpoint: pd.DataFrame,
    score_frame: pd.DataFrame,
    predictions: pd.DataFrame,
    model: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for score_name in score_frame.columns:
        scores = endpoint.loc[:, ["assessment_index", "y"]].copy()
        scores["cell_score"] = score_frame.loc[endpoint["assessment_index"], score_name].to_numpy(dtype=float)
        aligned = align_locked_predictions(scores, predictions, cohort, model)
        correlation = spearmanr(aligned["cell_score"], aligned["probability"])
        rows.append(
            {
                "cohort": cohort,
                "model": model,
                "method": "MCPcounter",
                "cell_score": score_name,
                "n": int(len(aligned)),
                "spearman_rho": float(correlation.statistic),
                "p_value_two_sided": float(correlation.pvalue),
                "mapping_status": "EVALUATED_EXACT_LABEL_AND_INDEX_MATCH",
                "descriptive_only": True,
            }
        )
    result = pd.DataFrame(rows)
    result["p_value_bh"] = multipletests(result["p_value_two_sided"], method="fdr_bh")[1]
    result["fdr_family"] = f"{cohort}|{model}|{len(result)}_scores"
    return result.sort_values("cell_score", kind="mergesort").reset_index(drop=True)


def add_cross_cohort_direction_metrics(summary: pd.DataFrame) -> pd.DataFrame:
    """Attach unpooled availability, direction counts, and an explicit concordance state."""
    valid_directions = {"MCI_higher", "HC_higher", "no_mean_difference"}
    if not summary["direction"].isin(valid_directions).all():
        raise ValueError("cross-cohort direction summary contains an invalid direction")
    rows: list[dict[str, object]] = []
    for cell_score, group in summary.groupby("cell_score", sort=True):
        counts = group["direction"].value_counts()
        available = int(len(group))
        mci_higher = int(counts.get("MCI_higher", 0))
        hc_higher = int(counts.get("HC_higher", 0))
        no_mean_difference = int(counts.get("no_mean_difference", 0))
        if mci_higher == available:
            concordance = "all_available_effects_MCI_higher"
        elif hc_higher == available:
            concordance = "all_available_effects_HC_higher"
        elif no_mean_difference == available:
            concordance = "all_available_effects_no_mean_difference"
        else:
            concordance = "mixed_direction_or_null"
        rows.append(
            {
                "cell_score": cell_score,
                "available_effect_count": available,
                "mci_higher_count": mci_higher,
                "hc_higher_count": hc_higher,
                "no_mean_difference_count": no_mean_difference,
                "direction_balance_mci_minus_hc": mci_higher - hc_higher,
                "direction_concordance": concordance,
            }
        )
    metrics = pd.DataFrame(rows)
    return summary.merge(metrics, on="cell_score", how="left", validate="many_to_one")


def main(output_dir: Path = OUT) -> int:
    all_scores, preassessment_manifest = _load_preassessment()
    source_registry = load_postlock_source_registry()
    validate_expression_sources(preassessment_manifest)
    validate_label_sources(source_registry)
    locked_predictions = load_locked_predictions(source_registry)
    summaries: list[pd.DataFrame] = []
    endpoint_cache: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    cohort_audit: dict[str, object] = {}
    for cohort in COHORTS:
        score_frame, labels = _cohort_scores_with_labels(all_scores, cohort)
        summaries.append(summarize_mci_hc_scores(cohort, score_frame, labels))
        endpoint = _endpoint_mapping(labels)
        endpoint_cache[cohort] = (endpoint, score_frame)
        cohort_audit[cohort] = {
            "all_locked_score_rows": int(len(score_frame)),
            "endpoint_rows": int(len(endpoint)),
            "label_counts_endpoint": {
                "HC": int((endpoint["y"] == 0).sum()),
                "MCI": int((endpoint["y"] == 1).sum()),
            },
        }

    summary = pd.concat(summaries, ignore_index=True)
    summary = add_cross_cohort_direction_metrics(summary)
    summary["cross_cohort_direction_summary"] = "descriptive_only_not_pooled"

    prediction_correlations: list[pd.DataFrame] = []
    prediction_audit: list[dict[str, object]] = []
    for cohort in COHORTS:
        endpoint, score_frame = endpoint_cache[cohort]
        for model in ["cross_study_stable_ridge", "pooled_ridge_nested_selector"]:
            try:
                prediction_correlations.append(
                    _describe_prediction_correlations(cohort, endpoint, score_frame, locked_predictions, model)
                )
                prediction_audit.append(
                    {"cohort": cohort, "model": model, "status": "EVALUATED_EXACT_LABEL_AND_INDEX_MATCH"}
                )
            except ValueError as error:
                prediction_audit.append(
                    {
                        "cohort": cohort,
                        "model": model,
                        "status": "NOT_EVALUATED_MAPPING_MISMATCH",
                        "reason": str(error),
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / SUMMARY_NAME, index=False)
    correlations = pd.concat(prediction_correlations, ignore_index=True) if prediction_correlations else pd.DataFrame()
    correlations.to_csv(output_dir / CORRELATIONS_NAME, index=False)
    audit = {
        "status": "POSTLOCK_SINGLE_METHOD_EXPLORATORY_COMPOSITION_AUDIT",
        "preassessment_scores_path": PREASSESSMENT_SCORES.name,
        "preassessment_scores_sha256": sha256_file(PREASSESSMENT_SCORES),
        "preassessment_manifest_sha256": sha256_file(PREASSESSMENT_MANIFEST),
        "postlock_source_registry_sha256": sha256_file(POSTLOCK_SOURCE_REGISTRY),
        "preassessment_diagnosis_labels_accessed": preassessment_manifest["diagnosis_labels_accessed"],
        "method": "MCPcounter",
        "xcell_status": preassessment_manifest["xcell"]["status"],
        "exploratory_only": True,
        "absolute_proportions_reported": False,
        "pooled_effect_reported": False,
        "diagnostic_model_modified": False,
        "statistics": {
            "score_effect": "MCI-minus-HC mean difference and Hedges_g",
            "group_test": "two-sided Welch t-test",
            "multiple_testing": "Benjamini-Hochberg within cohort and MCPcounter method across ten scores",
            "cross_cohort_summary": "available effects and explicit direction counts/concordance only; not pooled",
            "prediction_correlation": "descriptive Spearman correlation after exact label/index validation",
        },
        "cohort_audit": cohort_audit,
        "cross_cohort_direction_metrics": (
            summary.loc[
                :,
                [
                    "cell_score",
                    "available_effect_count",
                    "mci_higher_count",
                    "hc_higher_count",
                    "no_mean_difference_count",
                    "direction_balance_mci_minus_hc",
                    "direction_concordance",
                ],
            ]
            .drop_duplicates()
            .sort_values("cell_score", kind="mergesort")
            .to_dict(orient="records")
        ),
        "prediction_correlation_audit": prediction_audit,
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [
                PREASSESSMENT_SCORES,
                PREASSESSMENT_MANIFEST,
                POSTLOCK_SOURCE_REGISTRY,
                LOCKED_PREDICTIONS,
                LOCKED_PREDICTION_AUDIT,
                GSE63060_MATRIX,
                GSE60_LABELS,
                GSE63061_MATRIX,
                GSE61_LABELS,
                GSE140829_MATRIX,
                GPL15988_FULL,
                GSE249477_RAW,
                GSE249477_SERIES,
                PROTOCOL,
                Path(__file__).resolve(),
            ]
        },
    }
    (output_dir / AUDIT_NAME).write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": audit["status"],
                "summary_rows": int(len(summary)),
                "prediction_correlation_rows": int(len(correlations)),
                "prediction_evaluations": prediction_audit,
            }
        )
    )
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run post-lock exploratory MCPcounter composition assessment")
    parser.add_argument("--output-dir", type=Path, default=OUT)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
