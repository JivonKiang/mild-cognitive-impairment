from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from run_gse249477_rank_transport import load_target_expression, parse_assessment_labels
from run_pairwise_external import (
    GPL15988_FULL,
    GSE140829_MATRIX,
    GSE60_LABELS,
    GSE60_MATRIX,
    GSE61_LABELS,
    GSE61_MATRIX,
    bootstrap_auc,
    load_dev,
    load_external,
    load_external_expression_unlabeled,
    metric_row,
    sha256_file,
)
from transductive_alignment_core import assert_label_free_alignment_contract, coral_align_source_to_target


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
SEED = 20260814
LABELS = {"HC": 0, "MCI": 1}
CORAL_REGULARIZATION = 1e-3


def load_development() -> tuple[pd.DataFrame, np.ndarray, dict[str, int]]:
    x60, meta60 = load_dev(GSE60_MATRIX, GSE60_LABELS, "GSE63060")
    x61, meta61 = load_dev(GSE61_MATRIX, GSE61_LABELS, "GSE63061")
    common = sorted(set(x60.columns).intersection(x61.columns))
    expression = pd.concat([x60[common], x61[common]], axis=0)
    metadata = pd.concat([meta60, meta61], axis=0).loc[expression.index]
    keep = metadata["label"].isin(LABELS)
    retained = metadata.loc[keep]
    return (
        expression.loc[keep],
        retained["label"].map(LABELS).to_numpy(dtype=int),
        dict(sorted(retained["label"].value_counts().to_dict().items())),
    )


def predict_with_transductive_coral(
    source: pd.DataFrame, source_labels: np.ndarray, target: pd.DataFrame, pca_components: int
) -> tuple[np.ndarray, int]:
    """Fit source components, align latent covariances, then scale with target features only."""
    if list(source.columns) != list(target.columns):
        raise ValueError("source and target features must be identically ordered")
    component_count = min(pca_components, source.shape[0] - 1, target.shape[0] - 1, source.shape[1])
    if component_count < 2:
        raise ValueError("at least two PCA components are required for the CORAL audit")
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    source_imputed = imputer.fit_transform(source)
    target_imputed = imputer.transform(target)
    source_scaled = scaler.fit_transform(source_imputed)
    target_scaled = scaler.transform(target_imputed)
    projector = PCA(n_components=component_count, svd_solver="randomized", random_state=SEED)
    source_latent = projector.fit_transform(source_scaled)
    target_latent = projector.transform(target_scaled)
    aligned_source = coral_align_source_to_target(
        source_latent, target_latent, regularization=CORAL_REGULARIZATION
    )
    target_latent_scaler = StandardScaler().fit(target_latent)
    aligned_source = target_latent_scaler.transform(aligned_source)
    target_latent = target_latent_scaler.transform(target_latent)
    model = LogisticRegression(C=0.2, max_iter=1500, solver="lbfgs", random_state=SEED)
    model.fit(aligned_source, source_labels)
    return model.predict_proba(target_latent)[:, 1], int(component_count)


def metric_for_locked_assessment(
    y: np.ndarray, probability: np.ndarray, cohort: str, pca_components: int
) -> dict[str, Any]:
    row = metric_row(y, probability, "TRANSDUCTIVE_EXTERNAL_STRESS_TEST", "stage_mci_vs_hc", cohort, "transductive_coral_pca")
    row["auc_ci_low"], row["auc_ci_high"] = bootstrap_auc(y, probability, SEED + len(y))
    row["pca_components"] = int(pca_components)
    row["coral_regularization"] = CORAL_REGULARIZATION
    row["target_labels_used_for_alignment"] = False
    row["target_labels_used_for_model_selection"] = False
    row["target_labels_used_only_for_assessment"] = True
    return row


def main(output_dir: Path = OUT, pca_components: int = 50) -> int:
    development, source_labels, source_counts = load_development()
    target140, target140_audit = load_external_expression_unlabeled(
        GSE140829_MATRIX, GPL15988_FULL, set(development.columns)
    )
    target249, target249_ids, target249_audit = load_target_expression()
    shared_genes = sorted(set(development.columns).intersection(target140.columns).intersection(target249.columns))
    if len(shared_genes) != 11575:
        raise ValueError(f"expected 11,575 shared exact symbols, observed {len(shared_genes)}")

    source = development[shared_genes]
    probability140_all, components140 = predict_with_transductive_coral(
        source, source_labels, target140[shared_genes], pca_components
    )
    probability249_all, components249 = predict_with_transductive_coral(
        source, source_labels, target249[shared_genes], pca_components
    )
    contract = assert_label_free_alignment_contract(
        target_labels_used_for_alignment=False,
        target_labels_used_for_model_selection=False,
        target_labels_used_for_assessment=True,
    )

    # Prediction is locked above. Diagnosis metadata are loaded only for this final assessment block.
    _assessment_expression, meta140, _assessment_audit = load_external(
        GSE140829_MATRIX, GPL15988_FULL, set(development.columns)
    )
    mask140 = meta140["label"].isin(LABELS)
    y140 = meta140.loc[mask140, "label"].map(LABELS).to_numpy(dtype=int)
    probability140 = pd.Series(probability140_all, index=target140.index).loc[meta140.index[mask140]].to_numpy()

    labels249 = parse_assessment_labels(target249_ids)
    mask249 = labels249.isin(LABELS)
    y249 = labels249.loc[mask249].map(LABELS).to_numpy(dtype=int)
    probability249 = probability249_all[mask249.to_numpy()]

    metrics = [
        metric_for_locked_assessment(y140, probability140, "GSE140829", components140),
        metric_for_locked_assessment(y249, probability249, "GSE249477", components249),
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metrics).to_csv(output_dir / "transductive_coral_metrics.csv", index=False)
    pd.concat(
        [
            pd.DataFrame({"cohort": "GSE140829", "assessment_index": np.arange(len(y140)), "y": y140, "probability": probability140}),
            pd.DataFrame({"cohort": "GSE249477", "assessment_index": np.arange(len(y249)), "y": y249, "probability": probability249}),
        ],
        ignore_index=True,
    ).to_csv(output_dir / "transductive_coral_predictions_anonymized.csv", index=False)
    audit = {
        "status": "TRANSDUCTIVE_CORAL_RESULTS_GENERATED",
        "method": "source-fitted median imputation, StandardScaler and PCA followed by target-unlabeled latent CORAL alignment",
        "method_source": "Sun et al., Return of Frustratingly Easy Domain Adaptation, AAAI 2016, DOI 10.1609/aaai.v30i1.10306",
        "combat_status": "EXCLUDED_FROM_LABEL_FREE_ANALYSIS: target diagnosis is unavailable, so ComBat cannot protect the disease covariate during batch adjustment",
        "source_development_rows": int(len(source_labels)),
        "source_label_counts": source_counts,
        "shared_exact_symbol_count": int(len(shared_genes)),
        "pca_components_requested": int(pca_components),
        "coral_regularization": CORAL_REGULARIZATION,
        "post_coral_target_standardization": "target-unlabeled latent StandardScaler",
        "target_preassessment_rows": {"GSE140829": int(len(target140)), "GSE249477": int(len(target249))},
        "gse140829_unlabeled_input_audit": target140_audit,
        "gse249477_unlabeled_input_audit": target249_audit,
        "target_labels_used_for_alignment": False,
        "target_labels_used_for_model_selection": False,
        "target_labels_used_only_for_assessment": True,
        **contract,
        "allowed_claims": [
            "A batch-level transductive CORAL stress test used only target expression covariates before assessment.",
            "Target diagnosis labels were accessed only after predictions were locked for descriptive assessment.",
        ],
        "forbidden_claims": [
            "CORAL provides a deployable single-sample clinical model.",
            "CORAL establishes external validation, causal biology, or a generalizable MCI biomarker.",
            "ComBat was validly applied without target diagnosis information.",
        ],
        "source_hashes": {
            str(path): sha256_file(path)
            for path in [
                GSE60_MATRIX,
                GSE60_LABELS,
                GSE61_MATRIX,
                GSE61_LABELS,
                GSE140829_MATRIX,
                GPL15988_FULL,
                Path(__file__).resolve(),
                ROOT / "02_ANALYSIS" / "transductive_alignment_core.py",
            ]
        },
    }
    (output_dir / "transductive_coral_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(pd.DataFrame(metrics).to_json(orient="records", indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the target-unlabeled transductive CORAL stress test.")
    parser.add_argument("--output-dir", type=Path, default=OUT)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
