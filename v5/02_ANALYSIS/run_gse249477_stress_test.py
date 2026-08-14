from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "03_RESULTS"
RAW = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_raw_count_normalize_04-10-2025.csv.gz"
)
SERIES = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_series_matrix.txt.gz"
)
PANEL = [
    ("GPI", 1),
    ("PHF15", 1),
    ("VAT1", 1),
    ("EIF3B", 1),
    ("UCP2", 1),
    ("LAT1-3TM", 1),
    ("C22ORF39", -1),
    ("USP38", -1),
    ("ZSWIM6", -1),
    ("PRDM10", -1),
    ("RP2", -1),
    ("GCA", -1),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical(value: object) -> str:
    return str(value).strip().upper()


def parse_titles() -> pd.DataFrame:
    titles: list[str] = []
    with gzip.open(SERIES, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))[1:]
                break
    rows = []
    for title in titles:
        sample_match = re.search(r"\[(DK22119_[0-9]+)", title)
        if not sample_match:
            continue
        lower = title.casefold()
        if "mild cognitive impairment" in lower:
            label = "MCI"
        elif "cognitively normal" in lower:
            label = "HC"
        elif "alzheimer" in lower:
            label = "AD"
        else:
            label = "OTHER"
        age_match = re.search(r"([0-9]+)y\.o\.", lower)
        rows.append(
            {
                "sample": sample_match.group(1),
                "label": label,
                "age": float(age_match.group(1)) if age_match else np.nan,
            }
        )
    return pd.DataFrame(rows).drop_duplicates("sample").set_index("sample")


def bootstrap_auc(y: np.ndarray, score: np.ndarray, seed: int = 20260812) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    for _ in range(2000):
        index = rng.integers(0, len(y), len(y))
        if len(np.unique(y[index])) == 2:
            estimates.append(float(roc_auc_score(y[index], score[index])))
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def main() -> int:
    meta = parse_titles()
    raw = pd.read_csv(RAW, compression="gzip", low_memory=False)
    total_columns = [column for column in raw.columns if " - TPM" in column]
    sample_names = [re.match(r"(DK22119_[0-9]+)", column).group(1) for column in total_columns]
    values = raw.loc[:, total_columns].apply(pd.to_numeric, errors="coerce")
    values.columns = sample_names
    raw["symbol"] = raw["Name"].map(canonical)
    expression = {}
    for gene, _sign in PANEL:
        rows = values.loc[raw["symbol"] == gene]
        if rows.empty:
            continue
        # Zero TPM is a valid observed value; only absent exact symbols are missing.
        expression[gene] = np.log2(rows.add(1).mean(axis=0))
    expression = pd.DataFrame(expression).reindex(meta.index)
    present = [gene for gene, _sign in PANEL if gene in expression.columns]
    missing = [gene for gene, _sign in PANEL if gene not in expression.columns]
    # Standardize across all labelled samples before restricting the endpoint.
    # This prevents label-dependent endpoint filtering from changing preprocessing.
    standardized = (expression[present] - expression[present].mean()) / expression[present].std(ddof=0)
    keep = meta["label"].isin({"HC", "MCI"})
    meta = meta.loc[keep]
    standardized = standardized.loc[meta.index]
    score = pd.Series(0.0, index=meta.index)
    for gene, sign in PANEL:
        if gene in standardized:
            score = score + sign * standardized[gene]
    y = meta["label"].map({"HC": 0, "MCI": 1}).to_numpy(dtype=int)
    score_values = score.to_numpy(dtype=float)
    auc = float(roc_auc_score(y, score_values))
    ci_low, ci_high = bootstrap_auc(y, score_values)
    result = {
        "status": "PARTIAL_EXTERNAL_SENSITIVITY",
        "dataset": "GSE249477",
        "assay": "blood bulk RNA-seq",
        "role": "cross-platform stress test only",
        "sample_count_all_labels": int(len(parse_titles())),
        "hc_n": int((meta["label"] == "HC").sum()),
        "mci_n": int((meta["label"] == "MCI").sum()),
        "ad_n_excluded": int((parse_titles()["label"] == "AD").sum()),
        "present_locked_symbols": present,
        "missing_locked_symbols": missing,
        "exact_mapping": bool(not missing),
        "preprocessing": "log2(TPM+1), gene-wise standardization across all 62 labelled samples before endpoint restriction, locked signed score",
        "preprocessing_all_label_samples": int(len(parse_titles())),
        "target_labels_used_for_preprocessing": False,
        "auc_mci_vs_hc": auc,
        "auc_ci_low": ci_low,
        "auc_ci_high": ci_high,
        "score_contrast_mean_mci_minus_hc": float(score[meta["label"] == "MCI"].mean() - score[meta["label"] == "HC"].mean()),
        "source_hashes": {
            str(RAW): sha256_file(RAW),
            str(SERIES): sha256_file(SERIES),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gse249477_stress_test.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
