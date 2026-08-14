from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(r"E:\20241004_MCI\20260805_Human_Genomics\03_REPRODUCIBILITY")
SCORES = SOURCE_ROOT / "GSE270454_locked12_scores.csv"
MAPPING = SOURCE_ROOT / "GSE270454_locked12_mapping_audit.csv"
CONTRASTS = SOURCE_ROOT / "GSE270454_locked12_contrasts.csv"
OUT = ROOT / "03_RESULTS"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    scores = pd.read_csv(SCORES)
    mapping = pd.read_csv(MAPPING)
    contrasts = pd.read_csv(CONTRASTS)
    subset = scores[scores["stage"].isin(["HC", "MCI"])].copy()
    y = (subset["stage"] == "MCI").astype(int)
    auc = float(roc_auc_score(y, subset["score"]))
    result = {
        "status": "BOUNDARY_AUDIT_ONLY",
        "dataset": "GSE270454",
        "assay": "whole-blood bulk RNA-seq",
        "n_rows": int(len(scores)),
        "mci_vs_hc_n": int(len(subset)),
        "label_counts": {key: int(value) for key, value in scores["stage"].value_counts().to_dict().items()},
        "frozen_panel_mapping": mapping.iloc[0].to_dict(),
        "mci_vs_hc_row_level_auc": auc,
        "adjusted_contrasts": contrasts.to_dict(orient="records"),
        "source_hashes": {str(path): sha256(path) for path in [SCORES, MAPPING, CONTRASTS]},
        "subject_independence": "not established for a new LODO model; use as an independent frozen-score negative boundary",
        "allowed_claim": "The frozen score did not reproduce a positive MCI-versus-HC direction in this independent whole-blood audit.",
        "forbidden_claim": "This is not a validated complete-panel replication or a clinical validation study.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gse270454_boundary_audit.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
