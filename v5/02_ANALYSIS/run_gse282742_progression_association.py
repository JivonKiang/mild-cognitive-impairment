"""Exploratory, subject-level GSE282742 progression association analysis.

This module intentionally contains association and paired-change statistics only.
It neither trains nor assesses any predictive model and cannot provide causal evidence.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(r"E:\20241004_MCI\20260709\external_public_data_20260721_151749")
COUNTS_PATH = DATA_DIR / "GSE282742_Expected_count.full.gz"
METADATA_PATH = DATA_DIR / "GSE282742_sample_subject_audit.csv"
DEFAULT_OUTPUT_DIR = ROOT / "03_RESULTS"
REFERENCE_DIR = ROOT / "01_DATA_AUDIT" / "reference"
REACTOME_GMT_PATH = REFERENCE_DIR / "ReactomePathways_gmt_current" / "ReactomePathways.gmt"
DUAL_SOURCE_MAP_PATH = REFERENCE_DIR / "gse282742_dual_source_ensembl_symbol_map_20260813.csv"
GENE_SET_REGISTRY_PATH = ROOT / "00_PROTOCOL" / "gse282742_progression_gene_set_registry_v1.json"
PATHWAY_PERMUTATIONS = 1_000
PATHWAY_SEED = 28_274_213
PATHWAY_MIN_SIZE = 15
PATHWAY_MAX_SIZE = 500
PATHWAY_BOUNDARY = (
    "exploratory competitive pathway association from adjusted P-MCI-minus-S-MCI statistics; "
    "not independent predictive validation; not causal evidence"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_source_hashes(expected_hashes: dict[str, str]) -> None:
    for raw_path, expected_hash in expected_hashes.items():
        observed_hash = sha256_file(Path(raw_path))
        if observed_hash != expected_hash:
            raise ValueError(f"source SHA256 mismatch: {raw_path}")


def read_subject_metadata() -> pd.DataFrame:
    metadata = pd.read_csv(METADATA_PATH, dtype={"subject_id": str, "vgh_id": str})
    required = {"sample", "sex", "group", "subject_id", "age", "vgh_id"}
    if not required.issubset(metadata.columns):
        raise ValueError("GSE282742 metadata lacks required subject fields")
    if len(metadata) != 116 or metadata["vgh_id"].nunique() != 116:
        raise ValueError("GSE282742 metadata row or VGH identifier count changed")
    if metadata["vgh_id"].duplicated().any():
        raise ValueError("GSE282742 metadata contains duplicate VGH identifiers")
    return metadata.loc[:, ["sample", "sex", "group", "subject_id", "age", "vgh_id"]].copy()


def resolve_baseline_subjects(metadata: pd.DataFrame, enforce_cohort_counts: bool = True) -> pd.DataFrame:
    baseline = metadata.loc[metadata["group"].isin(["P-MCI", "S-MCI"])].copy()
    if enforce_cohort_counts and len(baseline) != 67:
        raise ValueError("GSE282742 baseline MCI row count changed")
    if baseline["subject_id"].isna().any():
        raise ValueError("GSE282742 baseline MCI rows lack subject identifiers")
    # GEO metadata provides no visit date. Preserve its original row order and retain
    # exactly one representative P-MCI/S-MCI record per subject without implying baseline.
    baseline = baseline.drop_duplicates("subject_id", keep="first").reset_index(drop=True)
    if enforce_cohort_counts and (len(baseline) != 52 or baseline["subject_id"].nunique() != 52):
        raise ValueError("GSE282742 does not resolve to 52 unique baseline subjects")
    if enforce_cohort_counts and baseline["group"].value_counts().to_dict() != {"P-MCI": 22, "S-MCI": 30}:
        raise ValueError("GSE282742 baseline group counts changed")
    return baseline


def complete_case_baseline(baseline: pd.DataFrame) -> pd.DataFrame:
    adjusted = baseline.dropna(subset=["age", "sex"]).copy()
    if len(adjusted) != 51:
        raise ValueError("GSE282742 adjusted baseline cohort must contain 51 complete cases")
    if adjusted["sex"].isin(["F", "M"]).all() is False:
        raise ValueError("GSE282742 sex coding changed")
    return adjusted.reset_index(drop=True)


def resolve_paired_transitions(metadata: pd.DataFrame, enforce_cohort_counts: bool = True) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for subject, subject_rows in metadata.groupby("subject_id", sort=True):
        p_mci = subject_rows.loc[subject_rows["group"] == "P-MCI"].sort_values("vgh_id", kind="stable")
        ad = subject_rows.loc[subject_rows["group"] == "AD"].sort_values("vgh_id", kind="stable")
        if len(p_mci) >= 1 and len(ad) >= 1:
            rows.append(
                {
                    "subject_id": str(subject),
                    "baseline_group": "P-MCI",
                    "transition_group": "AD",
                    "baseline_vgh_id": str(p_mci.iloc[0]["vgh_id"]),
                    "transition_vgh_id": str(ad.iloc[0]["vgh_id"]),
                }
            )
    pairs = pd.DataFrame(rows).sort_values("subject_id", kind="stable").reset_index(drop=True)
    if enforce_cohort_counts and (len(pairs) != 11 or pairs["subject_id"].nunique() != 11):
        raise ValueError("GSE282742 must contain exactly 11 P-MCI-to-AD pairs")
    return pairs


def read_counts() -> pd.DataFrame:
    counts = pd.read_csv(COUNTS_PATH, sep="\t", compression="gzip", dtype={"gene_id": str})
    if counts.shape[1] != 117 or "gene_id" not in counts.columns:
        raise ValueError("GSE282742 count matrix dimensions changed")
    counts = counts.set_index("gene_id")
    if counts.index.duplicated().any() or (counts < 0).any().any():
        raise ValueError("GSE282742 count matrix has duplicate genes or negative counts")
    return counts


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("GSE282742 includes an empty library")
    return np.log2(counts.div(library_sizes, axis=1) * 1_000_000 + 1.0)


def parse_gmt(path: Path) -> dict[str, set[str]]:
    pathways: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name, reactome_id, *genes = fields
            gene_set = {gene for gene in genes if gene}
            if gene_set:
                pathways[f"{name} [{reactome_id}]"] = gene_set
    if not pathways:
        raise ValueError("Reactome GMT contains no gene sets")
    return pathways


def enrichment_score(ranked_symbols: np.ndarray, member_symbols: set[str]) -> float:
    hits = np.isin(ranked_symbols, list(member_symbols))
    hit_count = int(hits.sum())
    miss_count = len(ranked_symbols) - hit_count
    if hit_count == 0 or miss_count == 0:
        raise ValueError("gene set has no valid enrichment walk")
    walk = np.cumsum(np.where(hits, 1.0 / hit_count, -1.0 / miss_count))
    maximum = float(walk.max())
    minimum = float(walk.min())
    return maximum if abs(maximum) >= abs(minimum) else minimum


def preranked_gene_set_enrichment(
    ranking: pd.DataFrame,
    pathways: dict[str, set[str]],
    *,
    seed: int,
    min_size: int,
    max_size: int,
) -> pd.DataFrame:
    required = {"gene_symbol", "t_statistic"}
    if not required.issubset(ranking.columns):
        raise ValueError("ranking must contain gene_symbol and t_statistic")
    ranked = ranking.loc[:, ["gene_symbol", "t_statistic"]].dropna().copy()
    ranked = ranked.drop_duplicates("gene_symbol", keep="first").sort_values(
        ["t_statistic", "gene_symbol"], ascending=[False, True], kind="stable"
    )
    ranked_symbols = ranked["gene_symbol"].to_numpy(dtype=str)
    universe = set(ranked_symbols)
    rank_values = stats.rankdata(ranked["t_statistic"].to_numpy(dtype=float), method="average")
    rows: list[dict[str, object]] = []
    for pathway_name in sorted(pathways):
        overlap = pathways[pathway_name] & universe
        if not (min_size <= len(overlap) <= max_size):
            continue
        observed = enrichment_score(ranked_symbols, overlap)
        member_mask = np.isin(ranked_symbols, list(overlap))
        member_ranks = rank_values[member_mask]
        background_ranks = rank_values[~member_mask]
        rank_test = stats.mannwhitneyu(member_ranks, background_ranks, alternative="two-sided", method="asymptotic")
        rank_biserial = (2.0 * rank_test.statistic / (len(member_ranks) * len(background_ranks))) - 1.0
        rows.append(
            {
                "analysis": "exploratory_preranked_reactome_enrichment",
                "pathway_name": pathway_name,
                "overlap_gene_count": len(overlap),
                "enrichment_score": observed,
                "rank_biserial_enrichment": rank_biserial,
                "p_value_two_sided": rank_test.pvalue,
                "test_method": "two-sided asymptotic Mann-Whitney competitive rank test",
                "random_seed": seed,
                "evidence_boundary": PATHWAY_BOUNDARY,
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("no Reactome pathways pass the configured overlap-size range")
    result["fdr_bh"] = multipletests(result["p_value_two_sided"], method="fdr_bh")[1]
    return result.sort_values(["fdr_bh", "p_value_two_sided", "pathway_name"], kind="stable").reset_index(drop=True)


def load_pinned_pathway_inputs() -> tuple[dict[str, object], pd.DataFrame, dict[str, set[str]]]:
    registry = json.loads(GENE_SET_REGISTRY_PATH.read_text(encoding="utf-8"))
    expected_gmt_hash = registry["reactome_gmt"]["gmt_sha256"]
    expected_map_hash = registry["mapping"]["sha256"]
    if sha256_file(REACTOME_GMT_PATH) != expected_gmt_hash:
        raise ValueError("Reactome GMT SHA256 does not match the pinned registry")
    if sha256_file(DUAL_SOURCE_MAP_PATH) != expected_map_hash:
        raise ValueError("dual-source Ensembl-symbol map SHA256 does not match the pinned registry")
    mapping = pd.read_csv(DUAL_SOURCE_MAP_PATH, dtype=str)
    if list(mapping.columns) != ["ensembl_gene_id", "gene_symbol"] or mapping.duplicated("ensembl_gene_id").any():
        raise ValueError("dual-source Ensembl-symbol map has an invalid schema")
    if len(mapping) != int(registry["mapping"]["row_count"]):
        raise ValueError("dual-source Ensembl-symbol map row count does not match the pinned registry")
    return registry, mapping, parse_gmt(REACTOME_GMT_PATH)


def build_pathway_ranking(baseline_results: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    ranking = baseline_results.loc[:, ["ensembl_gene_id", "t_statistic"]].merge(
        mapping, on="ensembl_gene_id", how="inner", validate="one_to_one"
    )
    ranking = ranking.dropna(subset=["t_statistic", "gene_symbol"]).copy()
    ranking["abs_t_statistic"] = ranking["t_statistic"].abs()
    ranking = ranking.sort_values(["gene_symbol", "abs_t_statistic", "ensembl_gene_id"], ascending=[True, False, True], kind="stable")
    ranking = ranking.drop_duplicates("gene_symbol", keep="first")
    return ranking.loc[:, ["gene_symbol", "t_statistic"]].sort_values(["t_statistic", "gene_symbol"], ascending=[False, True], kind="stable")


def fit_baseline_associations(log_expression: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    sample_ids = baseline["vgh_id"].tolist()
    if set(sample_ids) - set(log_expression.columns):
        raise ValueError("baseline metadata VGH identifiers are absent from the count matrix")
    design = pd.DataFrame(
        {
            "p_mci": (baseline["group"] == "P-MCI").astype(float).to_numpy(),
            "age": baseline["age"].astype(float).to_numpy(),
            "male": (baseline["sex"] == "M").astype(float).to_numpy(),
        }
    )
    design = sm.add_constant(design, has_constant="add")
    values = log_expression.loc[:, sample_ids].T.to_numpy(dtype=float)
    coefficients = np.linalg.lstsq(design.to_numpy(), values, rcond=None)[0]
    residuals = values - design.to_numpy() @ coefficients
    degrees_freedom = len(baseline) - design.shape[1]
    residual_variance = (residuals**2).sum(axis=0) / degrees_freedom
    covariance_diagonal = np.diag(np.linalg.inv(design.to_numpy().T @ design.to_numpy()))
    standard_errors = np.sqrt(residual_variance * covariance_diagonal[1])
    estimate = coefficients[1]
    statistic = np.divide(estimate, standard_errors, out=np.full_like(estimate, np.nan), where=standard_errors > 0)
    p_values = 2.0 * stats.t.sf(np.abs(statistic), degrees_freedom)
    confidence = stats.t.ppf(0.975, degrees_freedom) * standard_errors
    fdr = np.full(len(p_values), np.nan)
    valid = np.isfinite(p_values)
    fdr[valid] = multipletests(p_values[valid], method="fdr_bh")[1]
    return pd.DataFrame(
        {
            "ensembl_gene_id": log_expression.index,
            "analysis": "exploratory_progression_association_subject_level_representative_record",
            "contrast": "P-MCI_minus_S-MCI_adjusted_for_age_and_sex",
            "estimate_log2_cpm": estimate,
            "ci95_low": estimate - confidence,
            "ci95_high": estimate + confidence,
            "standard_error": standard_errors,
            "t_statistic": statistic,
            "p_value_two_sided": p_values,
            "fdr_bh": fdr,
            "n_subjects": len(baseline),
            "df_residual": degrees_freedom,
            "evidence_boundary": "exploratory progression association; not independent predictive validation; not causal evidence",
        }
    )


def fit_paired_transitions(log_expression: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    baseline_ids = pairs["baseline_vgh_id"].tolist()
    transition_ids = pairs["transition_vgh_id"].tolist()
    if (set(baseline_ids) | set(transition_ids)) - set(log_expression.columns):
        raise ValueError("paired metadata VGH identifiers are absent from the count matrix")
    deltas = log_expression.loc[:, transition_ids].to_numpy() - log_expression.loc[:, baseline_ids].to_numpy()
    zero_variance = np.isclose(deltas.std(axis=1, ddof=1), 0.0)
    testable = ~zero_variance
    p_values = np.full(deltas.shape[0], np.nan)
    test_statistics = np.full(deltas.shape[0], np.nan)
    if testable.any():
        test = stats.ttest_1samp(deltas[testable], popmean=0.0, axis=1, nan_policy="raise")
        p_values[testable] = test.pvalue
        test_statistics[testable] = test.statistic
    mean_delta = deltas.mean(axis=1)
    standard_errors = deltas.std(axis=1, ddof=1) / np.sqrt(deltas.shape[1])
    confidence = stats.t.ppf(0.975, deltas.shape[1] - 1) * standard_errors
    fdr = np.full(deltas.shape[0], np.nan)
    valid = np.isfinite(p_values)
    fdr[valid] = multipletests(p_values[valid], method="fdr_bh")[1]
    return pd.DataFrame(
        {
            "ensembl_gene_id": log_expression.index,
            "analysis": "paired_transition_sensitivity",
            "contrast": "AD_minus_P-MCI_within_subject",
            "mean_delta_log2_cpm": mean_delta,
            "ci95_low": mean_delta - confidence,
            "ci95_high": mean_delta + confidence,
            "t_statistic": test_statistics,
            "p_value_two_sided": p_values,
            "fdr_bh": fdr,
            "n_pairs": len(pairs),
            "test_status": np.where(valid, "TESTABLE", "NOT_TESTABLE_ZERO_VARIANCE"),
            "evidence_boundary": "paired transition sensitivity; not independent predictive validation; not causal evidence",
        }
    )


def main(output_dir: Path = DEFAULT_OUTPUT_DIR) -> int:
    metadata = read_subject_metadata()
    baseline = complete_case_baseline(resolve_baseline_subjects(metadata))
    pairs = resolve_paired_transitions(metadata)
    counts = read_counts()
    log_expression = log2_cpm(counts)
    baseline_results = fit_baseline_associations(log_expression, baseline)
    paired_results = fit_paired_transitions(log_expression, pairs)
    pathway_registry, pathway_mapping, reactome_pathways = load_pinned_pathway_inputs()
    pathway_ranking = build_pathway_ranking(baseline_results, pathway_mapping)
    pathway_summary = preranked_gene_set_enrichment(
        pathway_ranking,
        reactome_pathways,
        seed=int(pathway_registry["statistics"]["seed"]),
        min_size=int(pathway_registry["statistics"]["pathway_overlap_size_range"][0]),
        max_size=int(pathway_registry["statistics"]["pathway_overlap_size_range"][1]),
    )
    audit = {
        "status": "EXPLORATORY_GSE282742_PROGRESSION_ASSOCIATION",
        "dataset": "GSE282742",
        "primary_analysis": "subject-level representative-record analysis",
        "primary_counts": {"unique_subjects": 52, "covariate_complete_subjects": 51, "p_mci": 22, "s_mci": 29},
        "representative_record_rule": "first eligible P-MCI/S-MCI row in public metadata order; no verified visit dates are available, so this is not claimed as baseline",
        "paired_analysis": "paired transition sensitivity",
        "paired_p_mci_to_ad_subjects": 11,
        "covariates": ["age", "sex"],
        "transformation": "log2(CPM + 1)",
        "multiplicity": "Benjamini-Hochberg across all tested genes per analysis",
        "prediction_model_fitting": False,
        "pathway_summary_status": "EXPLORATORY_PINNED_REACTOME_PRERANKED_ENRICHMENT",
        "pathway_summary_source": {
            "registry_path": str(GENE_SET_REGISTRY_PATH),
            "registry_sha256": sha256_file(GENE_SET_REGISTRY_PATH),
            "gmt_sha256": sha256_file(REACTOME_GMT_PATH),
            "mapping_sha256": sha256_file(DUAL_SOURCE_MAP_PATH),
            "ranked_unique_symbols": len(pathway_ranking),
            "tested_pathways": len(pathway_summary),
        },
        "claim_boundary": "exploratory progression association; subject-level representative-record analysis without verified visit dates; paired transition sensitivity; not independent predictive validation; not causal evidence",
        "source_hashes": {str(path): sha256_file(path) for path in [COUNTS_PATH, METADATA_PATH]},
        "official_sources": {
            "geo_record": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282742",
            "pubmed": "https://pubmed.ncbi.nlm.nih.gov/39754267/",
            "doi": "https://doi.org/10.1186/s13195-024-01651-0",
            "access_date": "2026-08-13",
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline_results.to_csv(output_dir / "gse282742_progression_baseline_association.csv", index=False)
    paired_results.to_csv(output_dir / "gse282742_progression_paired_transition.csv", index=False)
    pathway_summary.to_csv(output_dir / "gse282742_progression_pathway_summary.csv", index=False)
    (output_dir / "gse282742_progression_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
