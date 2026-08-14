"""Re-run the frozen 12-gene GSE270454 audit with the verified GEO join.

This script uses the locked panel and fixed signed score only.  It does not
select genes, tune weights, choose thresholds, or use target labels for any
adaptation.  The verified count-column/GEO mapping supplies age and sex.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(r"E:\20241004_MCI\20260805_Human_Genomics\03_REPRODUCIBILITY")
COUNTS = DATA_ROOT / "GSE270454_RNAseq-combined-counts-matrix.csv.gz"
PANEL = Path(r"E:\20241004_MCI\20260709\mci_development_optimized_panel_20260722_161919\development_locked_12_gene_panel.csv")
JOIN = ROOT / "01_DATA_AUDIT" / "gse270454_join_mapping_local.csv"
OUT = ROOT / "03_RESULTS"
SOURCE_OUT = ROOT / "04_FIGURES" / "source_data"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    mean = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan)
    return expr.sub(mean, axis=0).div(sd, axis=0).fillna(0.0)


def bootstrap_auc(y: np.ndarray, score: np.ndarray, seed: int = 270454, repeats: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(repeats):
        indices = rng.integers(0, len(y), len(y))
        if len(np.unique(y[indices])) < 2:
            continue
        values.append(float(roc_auc_score(y[indices], score[indices])))
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def contrasts(meta: pd.DataFrame) -> pd.DataFrame:
    design = pd.DataFrame(
        {
            "const": 1.0,
            "MCI": (meta["stage"] == "MCI").astype(float),
            "AD": (meta["stage"] == "AD").astype(float),
            "age": meta["age"].astype(float),
        }
    )
    sex = pd.get_dummies(meta["sex"].astype(str), prefix="sex", drop_first=True, dtype=float)
    design = pd.concat([design, sex.reset_index(drop=True)], axis=1)
    fit = sm.OLS(meta["score"].to_numpy(float), design.to_numpy(float)).fit()
    names = list(design.columns)
    rows: list[dict[str, object]] = []
    for label, weights in {
        "MCI_vs_HC": {"MCI": 1.0},
        "AD_vs_MCI": {"AD": 1.0, "MCI": -1.0},
        "AD_vs_HC": {"AD": 1.0},
    }.items():
        vector = np.zeros(len(names))
        for name, weight in weights.items():
            vector[names.index(name)] = weight
        estimate = float(vector @ fit.params)
        standard_error = float(np.sqrt(vector @ fit.cov_params() @ vector))
        p_value = float(np.asarray(fit.t_test(vector).pvalue).squeeze())
        rows.append(
            {
                "dataset": "GSE270454",
                "contrast": label,
                "estimate": estimate,
                "se": standard_error,
                "p_value": p_value,
                "n": len(meta),
                "covariates": "age,sex",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    mapping = pd.read_csv(JOIN)
    if len(mapping) != 45 or not mapping["join_status"].eq("unique").all():
        raise RuntimeError("Verified GSE270454 join is missing, non-unique, or incomplete")
    if not mapping["whole_blood"].all():
        raise RuntimeError("Verified join contains a non-whole-blood sample")

    counts = pd.read_csv(COUNTS, compression="gzip", index_col=0)
    counts.columns = [str(column).strip() for column in counts.columns]
    sample_columns = mapping["count_column"].tolist()
    expr = counts[sample_columns].apply(pd.to_numeric, errors="coerce")
    expr.index = expr.index.astype(str).str.strip().str.upper()
    expr = expr[(expr.index != "") & (expr.index != "NAN")]
    expr = expr.groupby(level=0, sort=False).median()

    panel = pd.read_csv(PANEL)
    panel["gene"] = panel["gene"].astype(str).str.strip().str.upper()
    panel = panel.drop_duplicates("gene")
    mapped = [gene for gene in panel["gene"] if gene in expr.index]
    missing = [gene for gene in panel["gene"] if gene not in expr.index]
    use = panel[panel["gene"].isin(mapped)].copy()
    standardized = zscore_rows(expr.loc[mapped, sample_columns])
    peak = use[use["trajectory_class"].eq("MCI_peak")]["gene"].tolist()
    trough = use[use["trajectory_class"].eq("MCI_trough")]["gene"].tolist()
    weights = pd.Series(0.0, index=mapped)
    weights.loc[peak] = 0.5 / len(peak)
    weights.loc[trough] = -0.5 / len(trough)
    score = standardized.mul(weights, axis=0).sum(axis=0)

    stage_by_group = {"ASM": "HC", "ASO": "HC", "MCI": "MCI", "AD": "AD"}
    meta = mapping.copy()
    meta["stage"] = meta["count_group"].map(stage_by_group)
    meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["sex"] = meta["sex"].astype(str).str.strip().str.lower()
    meta["score"] = meta["count_column"].map(score)
    meta = meta[["count_column", "geo_accession", "geo_index", "count_group", "stage", "age", "sex", "score"]]

    subset = meta[meta["stage"].isin(["HC", "MCI"])].copy()
    y = (subset["stage"] == "MCI").astype(int).to_numpy()
    score_values = subset["score"].to_numpy(float)
    auc = float(roc_auc_score(y, score_values))
    auc_ci_low, auc_ci_high = bootstrap_auc(y, score_values)
    contrast_table = contrasts(meta.dropna(subset=["age", "score", "sex"]))
    stage_summary = meta.groupby("stage", observed=False)["score"].agg(["count", "mean", "std", "median"]).reset_index()

    mapping_summary = {
        "dataset": "GSE270454",
        "panel_size": len(panel),
        "mapped": len(mapped),
        "missing": len(missing),
        "mapping_fraction": len(mapped) / len(panel),
        "missing_genes": ";".join(missing),
        "rna_samples": len(meta),
        "HC": int((meta["stage"] == "HC").sum()),
        "MCI": int((meta["stage"] == "MCI").sum()),
        "AD": int((meta["stage"] == "AD").sum()),
        "join_status": "verified_unique_45_of_45",
        "join_audit": str(ROOT / "03_RESULTS" / "gse270454_join_audit.json"),
    }
    audit = {
        "status": "BOUNDARY_AUDIT_ONLY",
        "dataset": "GSE270454",
        "assay": "whole-blood bulk RNA-seq",
        "join_status": "VERIFIED_UNIQUE_45_OF_45",
        "join_method": "count-column group and ordinal to descriptive GEO title group and terminal sample ordinal",
        "n_rows": int(len(meta)),
        "mci_vs_hc_n": int(len(subset)),
        "label_counts": {key: int(value) for key, value in meta["stage"].value_counts().to_dict().items()},
        "frozen_panel_mapping": mapping_summary,
        "mci_vs_hc_row_level_auc": auc,
        "mci_vs_hc_auc_ci_low": auc_ci_low,
        "mci_vs_hc_auc_ci_high": auc_ci_high,
        "adjusted_contrasts": contrast_table.to_dict(orient="records"),
        "source_hashes": {
            str(COUNTS): sha256(COUNTS),
            str(PANEL): sha256(PANEL),
            str(JOIN): sha256(JOIN),
        },
        "subject_independence": "official series describes 45 adults; this audit does not construct a new longitudinal or LODO model",
        "allowed_claim": "The verified independent whole-blood frozen-score audit did not reproduce a positive MCI-versus-HC direction.",
        "forbidden_claim": "This is not a validated complete-panel replication, clinical validation, or causal study.",
        "online_metadata_verification": "verified against NCBI Entrez GEO summary for GSE270454 and PubMed PMID 40779089 on 2026-08-13",
    }

    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.mkdir(parents=True, exist_ok=True)
    meta.to_csv(OUT / "gse270454_joined_locked12_scores_local.csv", index=False)
    contrast_table.to_csv(OUT / "gse270454_joined_locked12_contrasts.csv", index=False)
    stage_summary.to_csv(OUT / "gse270454_joined_locked12_stage_summary.csv", index=False)
    (OUT / "gse270454_joined_mapping_audit.json").write_text(json.dumps(mapping_summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (OUT / "gse270454_boundary_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (SOURCE_OUT / "gse270454_boundary_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    pd.DataFrame(
        [
            {"metric": "mci_vs_hc_row_level_auc", "value": auc},
            {"metric": "mci_vs_hc_n", "value": len(subset)},
            {"metric": "panel_mapped", "value": len(mapped)},
            {"metric": "panel_size", "value": len(panel)},
            {"metric": "join_unique_count", "value": 45},
            {"metric": "join_status", "value": "verified_unique_45_of_45"},
        ]
    ).to_csv(SOURCE_OUT / "gse270454_boundary_summary.csv", index=False)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
