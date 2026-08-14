"""Source-gated sex-adjusted GSE208623 methylation and transcriptome concordance audit."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "01_DATA_AUDIT" / "gse208623"
RESULTS = ROOT / "03_RESULTS"
BETA_PATH = INPUT_DIR / "GSE208623_Norm_beta.csv.gz"
ANNOTATION_PATH = INPUT_DIR / "GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
SERIES_HEADER_PATH = INPUT_DIR / "GSE208623_series_matrix_header.gz"
SOURCE_MANIFEST_PATH = INPUT_DIR / "gse208623_official_source_manifest.json"
BETA_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208623/suppl/GSE208623_Norm_beta.csv.gz"
ANNOTATION_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
CROSS_COHORT_EFFECTS_PATH = RESULTS / "cross_cohort_effects.csv"
PROTOCOL = ROOT / "00_PROTOCOL" / "analysis_protocol_gse208623_cross_omics_v1.md"
BETA_EXPECTED_BYTES = 363_136_110
ANNOTATION_EXPECTED_BYTES = 169_013_044
LOCKED_TRANSCRIPTOME_REFERENCE_SHA256 = "aadeab45740f751e77f9e0294f03ddb600861f892db4415557cfabe3d7adbb41"


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    finite = np.isfinite(p_values)
    if not finite.any():
        return adjusted
    values = p_values[finite]
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    n = len(ranked)
    corrected = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.minimum(corrected, 1.0)
    adjusted[finite] = restored
    return adjusted


def beta_to_m_values(beta: np.ndarray) -> np.ndarray:
    beta = np.asarray(beta, dtype=float)
    if not np.isfinite(beta).all() or (beta < 0.0).any() or (beta > 1.0).any():
        raise ValueError("beta values must be finite and within [0, 1]")
    epsilon = 1e-6
    return np.log2((beta + epsilon) / (1.0 - beta + epsilon))


def fit_cpg_association(beta: np.ndarray, endpoint_mci: np.ndarray, sex_male: np.ndarray) -> dict[str, float]:
    beta = np.asarray(beta, dtype=float)
    endpoint_mci = np.asarray(endpoint_mci, dtype=float)
    sex_male = np.asarray(sex_male, dtype=float)
    if not (len(beta) == len(endpoint_mci) == len(sex_male)):
        raise ValueError("beta, endpoint and sex arrays must be aligned")
    if set(np.unique(endpoint_mci)) - {0.0, 1.0} or set(np.unique(sex_male)) - {0.0, 1.0}:
        raise ValueError("endpoint and sex must be binary")
    m_values = beta_to_m_values(beta)
    design = np.column_stack([np.ones(len(beta)), endpoint_mci, sex_male])
    coefficients, _, rank, _ = np.linalg.lstsq(design, m_values, rcond=None)
    if rank != 3:
        raise ValueError("sex-adjusted CpG association requires a full-rank intercept/MCI/sex design")
    degrees_freedom = len(beta) - rank
    if degrees_freedom <= 0:
        raise ValueError("insufficient degrees of freedom for sex-adjusted CpG association")
    residuals = m_values - design @ coefficients
    sigma_squared = float(residuals @ residuals / degrees_freedom)
    covariance = sigma_squared * np.linalg.inv(design.T @ design)
    standard_error = float(np.sqrt(covariance[1, 1]))
    statistic = float(coefficients[1] / standard_error) if standard_error > 0 else np.nan
    p_value = float(2.0 * stats.t.sf(abs(statistic), degrees_freedom)) if np.isfinite(statistic) else np.nan
    critical = float(stats.t.ppf(0.975, degrees_freedom))
    return {
        "mci_minus_hc_m_value": float(coefficients[1]),
        "standard_error": standard_error,
        "ci_lower": float(coefficients[1] - critical * standard_error),
        "ci_upper": float(coefficients[1] + critical * standard_error),
        "t_statistic": statistic,
        "p_value": p_value,
        "mci_minus_hc_mean_beta": float(beta[endpoint_mci == 1].mean() - beta[endpoint_mci == 0].mean()),
    }


def align_beta_columns(beta_columns: list[str], metadata: pd.DataFrame, expected_gsms: list[str]) -> pd.DataFrame:
    """Require beta columns to contain every official GSM exactly once and in metadata order."""
    if not {"assessment_index", "endpoint", "sex"}.issubset(metadata.columns):
        raise ValueError("GSE208623 metadata lacks required assessment/endpoint/sex columns")
    sample_columns = [column for column in beta_columns if column.startswith("GSM")]
    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError("GSE208623 beta columns are not one-to-one: duplicated GSM")
    if len(expected_gsms) != len(set(expected_gsms)) or len(expected_gsms) != len(metadata):
        raise ValueError("GSE208623 expected GSM records are not one-to-one with metadata")
    if set(sample_columns) != set(expected_gsms):
        raise ValueError("GSE208623 beta sample columns are not one-to-one with official GSM records")
    aligned = metadata.copy()
    aligned["assessment_index"] = pd.to_numeric(aligned["assessment_index"], errors="raise")
    if aligned["assessment_index"].duplicated().any() or sorted(aligned["assessment_index"].tolist()) != list(range(len(aligned))):
        raise ValueError("GSE208623 metadata assessment_index must be a unique consecutive integer sequence")
    aligned = aligned.sort_values("assessment_index", kind="mergesort").reset_index(drop=True)
    aligned["gsm"] = expected_gsms
    if not set(aligned["endpoint"]).issubset({"MCI", "HC", "AD"}):
        raise ValueError("GSE208623 endpoint metadata contains an unexpected label")
    if not set(aligned["sex"]).issubset({"female", "male"}):
        raise ValueError("GSE208623 sex metadata contains an unexpected label")
    return aligned


def unambiguous_gene_mapping(annotation: pd.DataFrame) -> pd.Series:
    """Return exact probe-to-symbol mappings only for single official symbol annotations."""
    required = {"IlmnID", "UCSC_RefGene_Name"}
    if not required.issubset(annotation.columns):
        raise ValueError(f"EPIC annotation lacks required columns: {sorted(required - set(annotation.columns))}")
    rows_by_probe: dict[str, list[set[str]]] = {}
    for probe, raw_symbols in annotation[["IlmnID", "UCSC_RefGene_Name"]].itertuples(index=False):
        probe = str(probe).strip()
        tokens = {token.strip().upper() for token in str(raw_symbols).split(";") if token.strip() and str(raw_symbols).lower() != "nan"}
        if probe:
            rows_by_probe.setdefault(probe, []).append(tokens)
    mappings = {
        probe: next(iter(row_tokens[0]))
        for probe, row_tokens in rows_by_probe.items()
        if all(len(tokens) == 1 for tokens in row_tokens)
        and len({next(iter(tokens)) for tokens in row_tokens}) == 1
    }
    return pd.Series(mappings, dtype="object", name="gene")


def parse_series_gsms() -> list[str]:
    """Extract the exact official GSE208623 GSM order from the source-locked header."""
    fields: dict[str, list[str]] = {}
    with gzip.open(SERIES_HEADER_PATH, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                key, *values = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                fields[key] = values
    gsms = fields.get("!Sample_geo_accession", [])
    if len(gsms) != 60 or len(set(gsms)) != 60:
        raise ValueError("source-locked GSE208623 metadata header does not contain 60 unique GSM records")
    return gsms


def official_subject_number_to_gsm(header_path: Path) -> tuple[list[str], dict[str, str]]:
    """Recover the official sample-number-to-GSM mapping from paired series-header fields."""
    fields: dict[str, list[str]] = {}
    with gzip.open(header_path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
                key, *values = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                fields[key] = values
    titles = fields.get("!Sample_title", [])
    gsms = fields.get("!Sample_geo_accession", [])
    if len(titles) != 60 or len(gsms) != 60 or len(set(gsms)) != 60:
        raise ValueError("official GSE208623 series header lacks 60 aligned title/GSM records")
    mapping: dict[str, str] = {}
    for title, gsm in zip(titles, gsms):
        match = re.search(r"\[([^\[\]]+)\]\s*$", title)
        if match is None:
            raise ValueError("official GSE208623 sample title lacks a terminal subject number")
        subject_number = match.group(1).strip()
        if not subject_number or subject_number in mapping:
            raise ValueError("official GSE208623 subject numbers are not one-to-one")
        mapping[subject_number] = gsm
    return gsms, mapping


def remap_official_beta_columns(beta_matrix: pd.DataFrame, header_path: Path) -> tuple[list[str], pd.DataFrame]:
    """Convert the official number-labelled beta columns into source-locked GSM columns."""
    if beta_matrix.shape[1] < 2:
        raise ValueError("processed beta matrix lacks a probe column and sample columns")
    expected_gsms, subject_number_to_gsm = official_subject_number_to_gsm(header_path)
    probe_column = beta_matrix.columns[0]
    sample_columns = beta_matrix.columns[1:].tolist()
    if len(sample_columns) != len(set(sample_columns)) or set(sample_columns) != set(subject_number_to_gsm):
        raise ValueError("official GSE208623 beta subject-number columns are not one-to-one with the series header")
    remapped = beta_matrix.rename(columns={probe_column: "IlmnID", **subject_number_to_gsm}).copy()
    if remapped.columns.duplicated().any():
        raise ValueError("official GSE208623 beta column remapping introduced duplicate names")
    return expected_gsms, remapped


def read_epic_annotation(path: Path) -> pd.DataFrame:
    """Read the official EPIC CSV after its Illumina metadata preamble."""
    header_line: int | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("IlmnID,"):
                if header_line is not None:
                    raise ValueError("official EPIC annotation has multiple IlmnID header rows")
                header_line = line_number
    if header_line is None:
        raise ValueError("official EPIC annotation lacks an IlmnID CSV header after the Illumina preamble")
    annotation = pd.read_csv(path, compression="gzip", skiprows=header_line, low_memory=False)
    _require_columns(annotation, {"IlmnID", "UCSC_RefGene_Name"}, "EPIC annotation")
    return annotation


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_locked_transcriptome_reference() -> None:
    """Reject a changed AddNeuroMed effect reference before concordance analysis."""
    observed = sha256_file(CROSS_COHORT_EFFECTS_PATH)
    if observed != LOCKED_TRANSCRIPTOME_REFERENCE_SHA256:
        raise ValueError(
            "locked transcriptome reference SHA-256 mismatch: "
            f"expected {LOCKED_TRANSCRIPTOME_REFERENCE_SHA256}, observed {observed}"
        )


def _require_columns(table: pd.DataFrame, required: set[str], table_name: str) -> None:
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{table_name} lacks required columns: {missing}")


def build_locked_transcriptome_reference(effects: pd.DataFrame) -> pd.Series:
    """Average the two prespecified MCI-minus-HC Hedges' g sources per gene."""
    _require_columns(
        effects,
        {"cohort", "gene", "hedges_g_mci_minus_hc", "status"},
        "cross-cohort effects",
    )
    cohorts = ["GSE63060", "GSE63061"]
    selected = effects.loc[effects["cohort"].isin(cohorts)].copy()
    selected["gene"] = selected["gene"].astype(str).str.strip().str.upper()
    selected["hedges_g_mci_minus_hc"] = pd.to_numeric(selected["hedges_g_mci_minus_hc"], errors="coerce")
    if selected.empty or not selected["status"].eq("estimable").all():
        raise ValueError("locked transcriptome reference requires estimable GSE63060 and GSE63061 effects")
    if (selected["gene"] == "").any() or not np.isfinite(selected["hedges_g_mci_minus_hc"]).all():
        raise ValueError("locked transcriptome reference contains missing or non-finite gene effects")
    if selected.duplicated(["cohort", "gene"]).any():
        raise ValueError("locked transcriptome reference has duplicate cohort-gene effects")
    counts = selected.groupby("gene", sort=False)["cohort"].nunique()
    if not counts.eq(2).all():
        raise ValueError("locked transcriptome reference requires exactly one effect from each prespecified cohort per gene")
    reference = selected.groupby("gene", sort=True)["hedges_g_mci_minus_hc"].mean()
    reference.name = "transcriptome_hedges_g_mci_minus_hc"
    return reference


def aggregate_gene_directional_z(cpg_associations: pd.DataFrame, probe_to_gene: pd.Series) -> pd.DataFrame:
    """Summarize eligible CpGs with signed, equal-weight Stouffer Z statistics."""
    _require_columns(cpg_associations, {"cpg", "mci_minus_hc_m_value", "p_value"}, "CpG associations")
    mapping = probe_to_gene.rename("gene").rename_axis("cpg").reset_index()
    merged = cpg_associations.merge(mapping, on="cpg", how="inner", validate="one_to_one")
    effects = pd.to_numeric(merged["mci_minus_hc_m_value"], errors="coerce").to_numpy(dtype=float)
    p_values = pd.to_numeric(merged["p_value"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(effects) & np.isfinite(p_values) & (p_values >= 0.0) & (p_values <= 1.0) & (effects != 0.0)
    merged = merged.loc[valid].copy()
    if merged.empty:
        return pd.DataFrame(columns=["gene", "probe_count", "methylation_z", "methylation_direction"])
    p_values = np.clip(p_values[valid], np.finfo(float).tiny, 1.0)
    merged["signed_z"] = np.sign(effects[valid]) * stats.norm.isf(p_values / 2.0)
    grouped = merged.groupby("gene", sort=True)["signed_z"]
    result = grouped.agg(probe_count="count", methylation_z=lambda values: float(values.sum() / np.sqrt(len(values)))).reset_index()
    result["methylation_direction"] = np.select(
        [result["methylation_z"] > 0.0, result["methylation_z"] < 0.0],
        ["MCI_higher_methylation", "HC_higher_methylation"],
        default="no_net_direction",
    )
    return result


def summarize_gene_concordance(
    methylation_genes: pd.DataFrame,
    transcriptome_reference: pd.Series,
    permutations: int = 10_000,
    rng: np.random.Generator | None = None,
) -> tuple[pd.DataFrame, dict[str, float | int]]:
    """Compute prespecified same-sign and opposite-sign summaries plus a label-free sign null."""
    _require_columns(methylation_genes, {"gene", "probe_count", "methylation_z"}, "methylation gene summary")
    if not isinstance(permutations, int) or permutations < 1:
        raise ValueError("permutations must be a positive integer")
    reference = pd.Series(transcriptome_reference, name="transcriptome_hedges_g_mci_minus_hc").copy()
    reference.index = reference.index.astype(str).str.strip().str.upper()
    if reference.index.has_duplicates or not np.isfinite(reference.to_numpy(dtype=float)).all():
        raise ValueError("locked transcriptome reference must have unique finite genes")
    merged = methylation_genes.merge(reference.rename_axis("gene").reset_index(), on="gene", how="inner", validate="one_to_one")
    merged = merged.sort_values("gene", kind="mergesort").reset_index(drop=True)
    if len(merged) < 2:
        raise ValueError("gene concordance requires at least two mapped genes in the locked reference overlap")
    methylation = merged["methylation_z"].to_numpy(dtype=float)
    expression = merged["transcriptome_hedges_g_mci_minus_hc"].to_numpy(dtype=float)
    spearman = float(stats.spearmanr(methylation, expression).statistic)
    methylation_direction = np.sign(methylation)
    transcriptome_direction = np.sign(expression)
    sign_eligible = (methylation_direction != 0.0) & (transcriptome_direction != 0.0)
    eligible_count = int(sign_eligible.sum())
    if eligible_count < 2:
        raise ValueError("gene sign concordance requires at least two nonzero-direction overlap genes")
    same_sign = methylation_direction[sign_eligible] == transcriptome_direction[sign_eligible]
    observed = float(same_sign.mean())
    generator = rng if rng is not None else np.random.default_rng(20_862_3)
    null = np.empty(permutations, dtype=float)
    for index in range(permutations):
        null[index] = (methylation_direction[sign_eligible] == generator.permutation(transcriptome_direction[sign_eligible])).mean()
    null_mean = float(null.mean())
    summary: dict[str, float | int] = {
        "overlap_gene_count": int(len(merged)),
        "sign_eligible_gene_count": eligible_count,
        "zero_or_no_direction_gene_count": int(len(merged) - eligible_count),
        "spearman_correlation": spearman,
        "same_sign_fraction": observed,
        "opposite_sign_fraction": float(1.0 - observed),
        "permutations": permutations,
        "gene_sign_permutation_null_mean": null_mean,
        "gene_sign_permutation_p_value": float((1 + np.count_nonzero(np.abs(null - null_mean) >= abs(observed - null_mean))) / (permutations + 1)),
    }
    merged["methylation_direction"] = np.where(methylation_direction > 0.0, "MCI_higher_methylation", np.where(methylation_direction < 0.0, "HC_higher_methylation", "no_net_direction"))
    merged["transcriptome_direction"] = np.where(transcriptome_direction > 0.0, "MCI_higher_expression", np.where(transcriptome_direction < 0.0, "HC_higher_expression", "no_net_direction"))
    merged["sign_eligible"] = sign_eligible
    merged["same_sign"] = pd.array([bool(left == right) if eligible else pd.NA for left, right, eligible in zip(methylation_direction, transcriptome_direction, sign_eligible)], dtype="boolean")
    return merged, summary


def run_cross_omics_association(
    beta_matrix: pd.DataFrame,
    annotation: pd.DataFrame,
    metadata: pd.DataFrame,
    expected_gsms: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the frozen MCI-versus-HC CpG associations and official one-to-one gene aggregation."""
    _require_columns(beta_matrix, {"IlmnID"}, "processed beta matrix")
    if beta_matrix["IlmnID"].isna().any() or beta_matrix["IlmnID"].astype(str).str.strip().eq("").any():
        raise ValueError("processed beta matrix IlmnID values must be nonempty")
    if beta_matrix["IlmnID"].duplicated().any():
        raise ValueError("processed beta matrix IlmnID values must be unique")
    _require_columns(metadata, {"assessment_index", "endpoint", "sex"}, "GSE208623 metadata")
    assessment_index = pd.to_numeric(metadata["assessment_index"], errors="coerce")
    if assessment_index.isna().any() or not np.array_equal(np.sort(assessment_index.to_numpy()), np.arange(60)):
        raise ValueError("GSE208623 metadata assessment_index must be unique and exactly 0..59")
    if len(metadata) != 60:
        raise ValueError("GSE208623 metadata must contain exactly 60 rows")
    endpoint_counts = metadata["endpoint"].value_counts().to_dict()
    if endpoint_counts != {"MCI": 20, "HC": 20, "AD": 20}:
        raise ValueError("GSE208623 metadata must contain exactly 20 MCI, 20 HC, and 20 AD rows")
    aligned = align_beta_columns(beta_matrix.columns.tolist(), metadata, expected_gsms)
    analysis = aligned.loc[aligned["endpoint"].isin(["MCI", "HC"])].copy()
    if len(analysis) != 40 or analysis["endpoint"].value_counts().to_dict() != {"MCI": 20, "HC": 20}:
        raise ValueError("GSE208623 model population must be exactly 20 MCI and 20 HC samples")
    if int((aligned["endpoint"] == "AD").sum()) != 20:
        raise ValueError("GSE208623 must exclude exactly 20 AD samples before model construction")
    analysis_gsms = analysis["gsm"].tolist()
    beta_values = beta_matrix.set_index("IlmnID")[analysis_gsms].apply(pd.to_numeric, errors="coerce")
    endpoint = (analysis["endpoint"] == "MCI").astype(float).to_numpy()
    sex_male = (analysis["sex"] == "male").astype(float).to_numpy()
    records: list[dict[str, float | int | str]] = []
    for cpg, values in beta_values.iterrows():
        values_array = values.to_numpy(dtype=float)
        try:
            result = fit_cpg_association(values_array, endpoint, sex_male)
        except (ValueError, np.linalg.LinAlgError):
            result = {
                "mci_minus_hc_m_value": np.nan, "standard_error": np.nan, "ci_lower": np.nan,
                "ci_upper": np.nan, "t_statistic": np.nan, "p_value": np.nan, "mci_minus_hc_mean_beta": np.nan,
            }
        result.update({"cpg": str(cpg), "analysis_n": len(analysis), "n_ad_excluded": int((aligned["endpoint"] == "AD").sum())})
        records.append(result)
    cpg_associations = pd.DataFrame(records)
    cpg_associations["p_value_bh"] = benjamini_hochberg(cpg_associations["p_value"].to_numpy(dtype=float))
    genes = aggregate_gene_directional_z(cpg_associations, unambiguous_gene_mapping(annotation))
    return cpg_associations, genes


def _validate_manifest_source(manifest: dict[str, object], logical_name: str, path: Path, url: str, expected_bytes: int) -> None:
    sources = manifest.get("sources")
    if not isinstance(sources, dict) or not isinstance(sources.get(logical_name), dict):
        raise ValueError(f"GSE208623 source manifest lacks required logical source: {logical_name}")
    record = sources[logical_name]
    if record.get("logical_name") != logical_name:
        raise ValueError(f"GSE208623 source manifest logical name mismatch for {logical_name}")
    if record.get("url") != url:
        raise ValueError(f"GSE208623 source manifest URL mismatch for {logical_name}")
    if record.get("expected_bytes") != expected_bytes:
        raise ValueError(f"GSE208623 source manifest expected byte mismatch for {logical_name}")
    if record.get("sha256") != sha256_file(path):
        raise ValueError(f"GSE208623 source manifest SHA256 mismatch for {logical_name}")
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\n\r")
    if record.get("header") != header:
        raise ValueError(f"GSE208623 source manifest header mismatch for {logical_name}")


def require_sources() -> None:
    required = {BETA_PATH: BETA_EXPECTED_BYTES, ANNOTATION_PATH: ANNOTATION_EXPECTED_BYTES, SERIES_HEADER_PATH: None, SOURCE_MANIFEST_PATH: None, CROSS_COHORT_EFFECTS_PATH: None, PROTOCOL: None}
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("GSE208623 cross-omics source gate failed: " + "; ".join(missing))
    wrong_size = [f"{path.name}: expected {size}, observed {path.stat().st_size}" for path, size in required.items() if size is not None and path.stat().st_size != size]
    if wrong_size:
        raise ValueError("GSE208623 cross-omics byte-size gate failed: " + "; ".join(wrong_size))
    for path in (BETA_PATH, ANNOTATION_PATH, SERIES_HEADER_PATH):
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
                if not handle.readline():
                    raise ValueError("empty gzip header")
        except (OSError, UnicodeError, EOFError, ValueError) as error:
            raise ValueError(f"GSE208623 source gzip gate failed for {path.name}: {error}") from error
    try:
        manifest = json.loads(SOURCE_MANIFEST_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"GSE208623 source manifest is unreadable: {SOURCE_MANIFEST_PATH.name}") from error
    if not isinstance(manifest, dict) or manifest.get("status") != "COMPLETE":
        raise ValueError("GSE208623 source manifest is incomplete and cannot be used for analysis")
    _validate_manifest_source(manifest, "beta", BETA_PATH, BETA_URL, BETA_EXPECTED_BYTES)
    _validate_manifest_source(manifest, "annotation", ANNOTATION_PATH, ANNOTATION_URL, ANNOTATION_EXPECTED_BYTES)
    require_locked_transcriptome_reference()


def main(output_dir: Path = RESULTS) -> int:
    require_sources()
    from audit_gse208623_metadata import parse_metadata

    metadata = pd.DataFrame(parse_metadata(SERIES_HEADER_PATH))
    beta_matrix = pd.read_csv(BETA_PATH, compression="gzip", low_memory=False)
    expected_gsms, beta_matrix = remap_official_beta_columns(beta_matrix, SERIES_HEADER_PATH)
    annotation = read_epic_annotation(ANNOTATION_PATH)
    cpg_associations, methylation_genes = run_cross_omics_association(beta_matrix, annotation, metadata, expected_gsms)
    transcriptome_reference = build_locked_transcriptome_reference(pd.read_csv(CROSS_COHORT_EFFECTS_PATH))
    concordance, concordance_summary = summarize_gene_concordance(methylation_genes, transcriptome_reference)
    output_dir.mkdir(parents=True, exist_ok=True)
    cpg_path = output_dir / "gse208623_cpg_associations.csv"
    gene_path = output_dir / "gse208623_gene_concordance.csv"
    audit_path = output_dir / "gse208623_cross_omics_source_audit.json"
    cpg_associations.to_csv(cpg_path, index=False)
    concordance.to_csv(gene_path, index=False)
    audit = {
        "status": "GSE208623_CROSS_OMICS_ASSOCIATION_AUDIT",
        "input_hashes": {
            "gse208623_beta": sha256_file(BETA_PATH),
            "gpl21145_annotation": sha256_file(ANNOTATION_PATH),
            "gse208623_series_matrix_header": sha256_file(SERIES_HEADER_PATH),
            "gse208623_official_source_manifest": sha256_file(SOURCE_MANIFEST_PATH),
            "cross_cohort_effects": sha256_file(CROSS_COHORT_EFFECTS_PATH),
            "analysis_protocol_gse208623_cross_omics": sha256_file(PROTOCOL),
        },
        "output_hashes": {
            "gse208623_cpg_associations": sha256_file(cpg_path),
            "gse208623_gene_concordance": sha256_file(gene_path),
        },
        "official_gsm_count": len(expected_gsms),
        "analysis_sample_count": int(cpg_associations["analysis_n"].iloc[0]),
        "ad_excluded_before_modeling": int(cpg_associations["n_ad_excluded"].iloc[0]),
        "cpg_count": int(len(cpg_associations)),
        "unambiguous_gene_count": int(len(methylation_genes)),
        "transcriptome_reference_cohorts": ["GSE63060", "GSE63061"],
        "concordance": concordance_summary,
        "evidence_boundary": "Sex-adjusted MCI-minus-HC CpG association and cross-omics concordance only; not prediction, model validation, MR, causal inference, IDAT processing, or pathway analysis.",
    }
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
