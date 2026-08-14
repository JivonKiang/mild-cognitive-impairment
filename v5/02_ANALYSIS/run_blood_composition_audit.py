"""Label-blind MCPcounter source-bridge preassessment scoring."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
import urllib.request
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from run_pairwise_external import load_gpl_symbol_map


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
R_SCRIPT = Path(r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe")
MCP_COUNTER_BRIDGE_TIMEOUT_SECONDS = 120

GSE63060_MATRIX = Path(r"E:\20241004_MCI\20250522\GEO\data\GSE63060_data_matrix.csv")
GSE63061_MATRIX = Path(r"E:\20241004_MCI\20260709\data\GSE63061_data_matrix.csv")
GSE140829_MATRIX = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse140829_task8\GSE140829_series_matrix.fresh.txt.gz"
)
GPL15988_FULL = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse140829_task9\GPL15988_full.txt"
)
GSE249477_RAW = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_raw_count_normalize_04-10-2025.csv.gz"
)

MCP_COUNTER_COMMIT = "b6eac73e91c246fcff0bb1a5c68a816cd588fc48"
MCP_COUNTER_SOURCE_BLOB_SHA = "bd6d2cc533f7635698eb870e8a3331a0a36fae82"
MCP_COUNTER_GENES_BLOB_SHA = "a37f70367d3dc83260e966e651c5d618dbf84c7c"
MCP_COUNTER_SOURCE_API_URL = (
    "https://api.github.com/repos/ebecht/MCPcounter/git/blobs/"
    + MCP_COUNTER_SOURCE_BLOB_SHA
)
MCP_COUNTER_GENES_API_URL = (
    "https://api.github.com/repos/ebecht/MCPcounter/git/blobs/"
    + MCP_COUNTER_GENES_BLOB_SHA
)
GITHUB_HEADERS = {
    "Accept": "application/vnd.github+json",
    "User-Agent": "codex-mcpcounter-source-bridge",
}


def sha256_file(path: Path) -> str:
    """Return the SHA256 digest of one declared raw input."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_symbol(value: object) -> str:
    return str(value).strip().strip('"').upper()


def _numeric(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _require_existing(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("declared MCPcounter bridge inputs are missing: " + ", ".join(missing))


def decode_github_blob_payload(payload: Mapping[str, Any], expected_sha: str) -> bytes:
    """Decode only the expected base64 GitHub blob payload."""
    if payload.get("sha") != expected_sha:
        raise ValueError("GitHub blob SHA differs from the pinned expected SHA")
    if payload.get("encoding") != "base64":
        raise ValueError("GitHub blob must use base64 encoding")
    content = payload.get("content")
    if not isinstance(content, str):
        raise ValueError("GitHub base64 blob content is absent")
    try:
        return base64.b64decode("".join(content.split()).encode("ascii"), validate=True)
    except (UnicodeEncodeError, ValueError) as error:
        raise ValueError("GitHub blob has invalid base64 content") from error


def fetch_pinned_github_blob(api_url: str, expected_sha: str) -> bytes:
    """Fetch one pinned source artifact through the GitHub blob API."""
    request = urllib.request.Request(api_url, headers=GITHUB_HEADERS)
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.load(response)
    if not isinstance(payload, dict):
        raise ValueError("GitHub API response must be a JSON object")
    return decode_github_blob_payload(payload, expected_sha)


def assert_expected_expression_rows(
    expression: pd.DataFrame,
    expected_count: int,
    cohort: str,
) -> None:
    """Require a usable, complete all-row expression frame before scoring."""
    if not isinstance(expression, pd.DataFrame):
        raise TypeError(f"{cohort} expression must be a pandas DataFrame")
    if expression.index.has_duplicates:
        raise ValueError(f"{cohort} expression row IDs must be unique")
    if len(expression) != expected_count:
        raise ValueError(f"{cohort} expected {expected_count} rows, found {len(expression)}")
    if expression.shape[1] == 0:
        raise ValueError(f"{cohort} expression contains no genes")
    if any(not str(gene).strip() for gene in expression.columns):
        raise ValueError(f"{cohort} expression contains a blank gene identifier")
    values = expression.to_numpy(dtype=float)
    if not np.isfinite(values).any():
        raise ValueError(f"{cohort} expression is all non-finite")


def _mean_duplicate_gene_rows(
    sums: dict[str, np.ndarray],
    counts: dict[str, np.ndarray],
    sample_ids: list[str],
) -> pd.DataFrame:
    if not sums:
        return pd.DataFrame(index=sample_ids)
    values = {
        gene: np.divide(
            sums[gene],
            counts[gene],
            out=np.full(len(sample_ids), np.nan, dtype=np.float64),
            where=counts[gene] > 0,
        )
        for gene in sorted(sums)
    }
    return pd.DataFrame(values, index=sample_ids)


def accumulate_expression_gene_row(
    raw_gene: object,
    values: list[float] | np.ndarray,
    sums: dict[str, np.ndarray],
    counts: dict[str, np.ndarray],
    cohort: str,
) -> None:
    """Accumulate one raw gene row while rejecting an ambiguous identifier."""
    gene = _canonical_symbol(raw_gene)
    if not gene:
        raise ValueError(f"{cohort} expression contains a blank gene identifier")
    array = np.asarray(values, dtype=np.float64)
    valid = np.isfinite(array)
    if gene not in sums:
        sums[gene] = np.zeros(len(array), dtype=np.float64)
        counts[gene] = np.zeros(len(array), dtype=np.float64)
    sums[gene][valid] += array[valid]
    counts[gene][valid] += 1.0


def load_development_expression_label_blind(matrix_path: Path, expected_count: int, cohort: str) -> pd.DataFrame:
    """Read every development matrix row without metadata or diagnosis access."""
    with matrix_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"{cohort} matrix is empty") from error
        if not header or header[0] != "Name":
            raise ValueError(f"{cohort} must have Name as its first column")
        sample_ids = [str(value).strip().strip('"') for value in header[1:]]
        if not sample_ids or len(set(sample_ids)) != len(sample_ids):
            raise ValueError(f"{cohort} expression row IDs must be nonempty and unique")

        sums: dict[str, np.ndarray] = {}
        counts: dict[str, np.ndarray] = {}
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"{cohort} matrix row has an unexpected field count")
            values = np.asarray([_numeric(value) for value in row[1:]], dtype=np.float64)
            accumulate_expression_gene_row(row[0], values, sums, counts, cohort)

    expression = _mean_duplicate_gene_rows(sums, counts, sample_ids)
    assert_expected_expression_rows(expression, expected_count, cohort)
    return expression


def load_gse140829_expression_label_blind() -> pd.DataFrame:
    """Map all valid platform HUGO symbols while reading only the matrix table."""
    probe_symbols = load_gpl_symbol_map(GPL15988_FULL)
    with gzip.open(GSE140829_MATRIX, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        sample_ids: list[str] = []
        sums: dict[str, np.ndarray] = {}
        counts: dict[str, np.ndarray] = {}
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line != "!series_matrix_table_begin":
                continue
            header = next(csv.reader([next(handle).rstrip("\r\n")], delimiter="\t"))
            sample_ids = [str(value).strip().strip('"') for value in header[1:]]
            if not sample_ids or len(set(sample_ids)) != len(sample_ids):
                raise ValueError("GSE140829 expression row IDs must be nonempty and unique")
            break
        if not sample_ids:
            raise ValueError("GSE140829 matrix table header was not found")

        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_end":
                break
            fields = next(csv.reader([line], delimiter="\t"))
            if len(fields) != len(sample_ids) + 1:
                raise ValueError("GSE140829 matrix row has an unexpected field count")
            symbols = probe_symbols.get(str(fields[0]).strip().strip('"'), set())
            if not symbols:
                continue
            values = np.asarray([_numeric(value) for value in fields[1:]], dtype=np.float64)
            valid = np.isfinite(values)
            for gene in symbols:
                if gene not in sums:
                    sums[gene] = np.zeros(len(sample_ids), dtype=np.float64)
                    counts[gene] = np.zeros(len(sample_ids), dtype=np.float64)
                sums[gene][valid] += values[valid]
                counts[gene][valid] += 1.0

    expression = _mean_duplicate_gene_rows(sums, counts, sample_ids)
    assert_expected_expression_rows(expression, 587, "GSE140829")
    return expression


def load_gse249477_expression_label_blind() -> tuple[pd.DataFrame, dict[str, object]]:
    """Load GSE249477 TPM data without opening any phenotype metadata file."""
    raw = pd.read_csv(
        GSE249477_RAW,
        compression="gzip",
        low_memory=False,
        keep_default_na=False,
        na_filter=False,
    )
    if "Name" not in raw.columns:
        raise ValueError("GSE249477 raw TPM table lacks the Name column")
    tpm_columns = [column for column in raw.columns if " - TPM" in column]
    if len(tpm_columns) != 62:
        raise ValueError(f"expected 62 GSE249477 TPM columns, found {len(tpm_columns)}")
    sample_ids: list[str] = []
    for column in tpm_columns:
        match = re.fullmatch(r"(DK22119_\d+) \(GE\) - TPM", column)
        if match is None:
            raise ValueError(f"GSE249477 TPM column has an unexpected header: {column!r}")
        sample_ids.append(match.group(1))
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("GSE249477 TPM identifiers are not unique")

    values = raw.loc[:, tpm_columns].apply(pd.to_numeric, errors="coerce")
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, np.ndarray] = {}
    for gene, row_values in zip(raw["Name"], values.to_numpy(dtype=np.float64), strict=True):
        accumulate_expression_gene_row(gene, row_values, sums, counts, "GSE249477")
    expression = _mean_duplicate_gene_rows(sums, counts, sample_ids)
    assert_expected_expression_rows(expression, 62, "GSE249477")
    audit = {
        "metadata_files_accessed": [],
        "target_tpm_columns": int(len(tpm_columns)),
        "unique_tpm_ids": int(len(set(sample_ids))),
        "target_gene_count_after_exact_symbol_mean_collapse": int(expression.shape[1]),
        "target_duplicate_gene_rows_collapsed": int(len(raw) - expression.shape[1]),
    }
    return expression.astype(np.float32), audit


def _read_transient_score_table(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("MCPcounter bridge did not write a score table") from error
        if len(header) < 2:
            raise ValueError("MCPcounter bridge score table has no score columns")
        row_ids: list[str] = []
        rows: list[list[float]] = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError("MCPcounter bridge score row has an unexpected field count")
            row_ids.append(row[0])
            rows.append([_numeric(value) for value in row[1:]])
    return pd.DataFrame(rows, index=row_ids, columns=header[1:])


def build_mcpcounter_bridge_script() -> str:
    """Build the disposable R bridge without changing the pinned upstream method."""
    return (
        "args <- commandArgs(trailingOnly=TRUE)\n"
        "source(args[[1]])\n"
        "genes <- read.table(args[[2]], sep='\\t', header=TRUE, stringsAsFactors=FALSE, "
        "colClasses='character', check.names=FALSE, quote='\\\"')\n"
        "expression <- read.table(args[[3]], sep='\\t', header=TRUE, row.names=1, "
        "check.names=FALSE, quote='', comment.char='', stringsAsFactors=FALSE, na.strings=character())\n"
        "markers <- split(as.character(genes[['HUGO symbols']]), "
        "as.character(genes[['Cell population']]))\n"
        "scores <- appendSignatures(expression, markers)\n"
        "write.table(scores, file=args[[4]], sep='\\t', quote=FALSE, col.names=NA)\n"
    )


def score_with_transient_mcpcounter_bridge(
    expression: pd.DataFrame,
    source_blob: bytes,
    genes_blob: bytes,
) -> pd.DataFrame:
    """Call the pinned upstream appendSignatures function in a disposable directory."""
    with tempfile.TemporaryDirectory(prefix="mcpcounter-source-bridge-") as directory:
        temporary = Path(directory)
        source_path = temporary / "MCPcounter.R"
        genes_path = temporary / "genes.txt"
        input_path = temporary / "expression.tsv"
        r_script_path = temporary / "run_append_signatures.R"
        scores_path = temporary / "scores.tsv"
        source_path.write_bytes(source_blob)
        genes_path.write_bytes(genes_blob)
        expression.T.to_csv(input_path, sep="\t", index=True, index_label="HUGO_symbol")
        r_script_path.write_text(build_mcpcounter_bridge_script(), encoding="utf-8")
        try:
            completed = subprocess.run(
                [str(R_SCRIPT), str(r_script_path), str(source_path), str(genes_path), str(input_path), str(scores_path)],
                capture_output=True,
                text=True,
                check=False,
                timeout=MCP_COUNTER_BRIDGE_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                "pinned MCPcounter appendSignatures bridge timed out after "
                f"{MCP_COUNTER_BRIDGE_TIMEOUT_SECONDS} seconds; "
                f"stdout={error.stdout!r}; stderr={error.stderr!r}"
            ) from error
        if completed.returncode != 0:
            raise RuntimeError(
                "pinned MCPcounter appendSignatures bridge failed with exit code "
                f"{completed.returncode}; stdout={completed.stdout[-1000:]!r}; stderr={completed.stderr[-2000:]!r}"
            )
        if not scores_path.exists():
            raise RuntimeError("pinned MCPcounter appendSignatures bridge produced no score table")
        return _read_transient_score_table(scores_path)


def score_all_rows_label_blind(
    expression: pd.DataFrame,
    scorer: Callable[[pd.DataFrame], pd.DataFrame],
) -> pd.DataFrame:
    """Score every expression row while preserving its exact row order."""
    if not isinstance(expression, pd.DataFrame):
        raise TypeError("expression must be a pandas DataFrame")
    if expression.index.has_duplicates:
        raise ValueError("expression rows must be unique")
    if not callable(scorer):
        raise TypeError("scorer must be callable")

    scores = scorer(expression.copy(deep=True))

    if not isinstance(scores, pd.DataFrame):
        raise TypeError("scorer must return a pandas DataFrame")
    if scores.index.has_duplicates:
        raise ValueError("scored rows must be unique")
    if scores.shape[1] == 0:
        raise ValueError("scorer must return at least one score column")
    if not scores.index.equals(expression.index):
        raise ValueError("scored rows must exactly match expression rows and order")
    if scores.isna().to_numpy().any():
        raise ValueError("scored values must not be missing")
    numeric_scores = scores.select_dtypes(include="number")
    if not numeric_scores.empty and not np.isfinite(numeric_scores.to_numpy()).all():
        raise ValueError("scored numeric values must be finite")

    return scores.loc[expression.index]


def assert_complete_preassessment_rows(
    expected_rows: Iterable[str],
    scored_rows: Iterable[str],
) -> None:
    """Require scored rows to be a duplicate-free match of expected rows."""
    expected = list(expected_rows)
    scored = list(scored_rows)

    if len(set(expected)) != len(expected):
        raise ValueError("expected preassessment rows contain duplicates")
    if len(set(scored)) != len(scored):
        raise ValueError("scored preassessment rows contain duplicates")

    expected_set = set(expected)
    scored_set = set(scored)
    missing = expected_set - scored_set
    extra = scored_set - expected_set
    if missing or extra:
        details = []
        if missing:
            details.append("missing: " + ", ".join(map(str, sorted(missing, key=str))))
        if extra:
            details.append("extra: " + ", ".join(map(str, sorted(extra, key=str))))
        raise ValueError("preassessment row membership differs (" + "; ".join(details) + ")")


def load_all_expression_label_blind() -> tuple[dict[str, pd.DataFrame], dict[str, object]]:
    """Load every declared expression row without reading phenotype metadata."""
    required = [
        GSE63060_MATRIX,
        GSE63061_MATRIX,
        GSE140829_MATRIX,
        GPL15988_FULL,
        GSE249477_RAW,
        R_SCRIPT,
    ]
    _require_existing(required)
    gse63060 = load_development_expression_label_blind(GSE63060_MATRIX, 329, "GSE63060")
    gse63061 = load_development_expression_label_blind(GSE63061_MATRIX, 388, "GSE63061")
    gse140829 = load_gse140829_expression_label_blind()
    gse249477, gse249477_audit = load_gse249477_expression_label_blind()
    expressions = {
        "GSE63060": gse63060,
        "GSE63061": gse63061,
        "GSE140829": gse140829,
        "GSE249477": gse249477,
    }
    provenance = {
        "input_hashes": {str(path): sha256_file(path) for path in required[:-1]},
        "expression_gene_counts": {cohort: int(frame.shape[1]) for cohort, frame in expressions.items()},
        "gse249477_label_free_loader_audit": gse249477_audit,
    }
    return expressions, provenance


def anonymize_preassessment_scores(
    cohort: str,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    """Drop source row IDs after enforcing complete preassessment coverage."""
    anonymized = scores.reset_index(drop=True).copy()
    anonymized.insert(0, "assessment_index", np.arange(len(anonymized), dtype=int))
    anonymized.insert(0, "cohort", cohort)
    return anonymized


def main(output_dir: Path = OUT) -> int:
    """Generate label-blind MCPcounter source-bridge scores for all cohort rows."""
    source_blob = fetch_pinned_github_blob(MCP_COUNTER_SOURCE_API_URL, MCP_COUNTER_SOURCE_BLOB_SHA)
    genes_blob = fetch_pinned_github_blob(MCP_COUNTER_GENES_API_URL, MCP_COUNTER_GENES_BLOB_SHA)
    expressions, provenance = load_all_expression_label_blind()

    expected_counts = {"GSE63060": 329, "GSE63061": 388, "GSE140829": 587, "GSE249477": 62}
    anonymized_parts: list[pd.DataFrame] = []
    cohort_audit: dict[str, object] = {}
    for cohort, expression in expressions.items():
        assert_expected_expression_rows(expression, expected_counts[cohort], cohort)
        scores = score_all_rows_label_blind(
            expression,
            lambda values: score_with_transient_mcpcounter_bridge(values, source_blob, genes_blob),
        )
        assert_complete_preassessment_rows(expression.index.tolist(), scores.index.tolist())
        cohort_audit[cohort] = {
            "expected_rows": expected_counts[cohort],
            "scored_rows": int(len(scores)),
            "score_column_count": int(scores.shape[1]),
        }
        anonymized_parts.append(anonymize_preassessment_scores(cohort, scores))

    anonymized = pd.concat(anonymized_parts, ignore_index=True)
    expected_total = sum(expected_counts.values())
    if len(anonymized) != expected_total:
        raise RuntimeError(f"expected {expected_total} anonymized score rows, found {len(anonymized)}")
    observed_counts = anonymized.groupby("cohort", sort=False).size().to_dict()
    if observed_counts != expected_counts:
        raise RuntimeError(f"unexpected anonymized cohort counts: {observed_counts}")
    if anonymized.columns[:2].tolist() != ["cohort", "assessment_index"]:
        raise RuntimeError("anonymized score table must begin with cohort and assessment_index")

    output_dir.mkdir(parents=True, exist_ok=True)
    scores_path = output_dir / "blood_composition_preassessment_scores_anonymized.csv"
    manifest_path = output_dir / "blood_composition_preassessment_manifest.json"
    anonymized.to_csv(scores_path, index=False)
    manifest = {
        "status": "MCP_COUNTER_SOURCE_BRIDGE_PREASSESSMENT",
        "method": "MCPcounter",
        "assessment_scope": "all expression rows before diagnosis assessment",
        "diagnosis_labels_accessed": False,
        "target_labels_used_for_mapping": False,
        "target_labels_used_for_scoring": False,
        "target_labels_used_for_preprocessing": False,
        "target_labels_used_for_method_selection": False,
        "diagnostic_model_modified": False,
        "xcell": {
            "status": "BLOCKED",
            "reason": "project-local xCell installation remains unavailable because R dependency retrieval is TLS-blocked",
        },
        "mcpcounter_source_bridge": {
            "repository": "https://github.com/ebecht/MCPcounter",
            "commit": MCP_COUNTER_COMMIT,
            "source_api_url": MCP_COUNTER_SOURCE_API_URL,
            "source_blob_sha": MCP_COUNTER_SOURCE_BLOB_SHA,
            "source_blob_sha256": hashlib.sha256(source_blob).hexdigest(),
            "genes_api_url": MCP_COUNTER_GENES_API_URL,
            "genes_blob_sha": MCP_COUNTER_GENES_BLOB_SHA,
            "genes_blob_sha256": hashlib.sha256(genes_blob).hexdigest(),
            "upstream_function_invoked": "appendSignatures",
            "bridge_implementation": "unmodified upstream function and fixed CC0 genes.txt sourced only in a disposable temporary directory",
            "upstream_source_or_signatures_persisted_in_project": False,
        },
        "input_hashes": provenance["input_hashes"],
        "expression_gene_counts": provenance["expression_gene_counts"],
        "gse249477_label_free_loader_audit": provenance["gse249477_label_free_loader_audit"],
        "cohort_preassessment": cohort_audit,
        "anonymized_output": {
            "path": scores_path.name,
            "sha256": sha256_file(scores_path),
            "row_count": int(len(anonymized)),
            "columns": anonymized.columns.tolist(),
            "original_sample_identifiers_persisted": False,
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "rows": int(len(anonymized)),
                "cohort_rows": observed_counts,
                "xcell_status": manifest["xcell"]["status"],
            }
        )
    )
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run label-blind MCPcounter source-bridge preassessment scoring")
    parser.add_argument("--output-dir", type=Path, default=OUT)
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main(parse_arguments().output_dir))
