"""Replay and audit the archived development-only locked 12-gene selection."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "03_RESULTS"
HISTORICAL_ROOT = Path(r"E:\20241004_MCI\20260709")
CANDIDATE_REGISTRY = HISTORICAL_ROOT / "mci_centered_score_denovo_20260722_140059" / "mci_centered_candidate_registry_weighted.csv"
BOOTSTRAP_STABILITY = HISTORICAL_ROOT / "mci_centered_panel_compression_20260722_144351" / "bootstrap_gene_stability.csv"
PROVISIONAL_PANELS = HISTORICAL_ROOT / "mci_centered_panel_compression_20260722_144351" / "provisional_panels_by_size.csv"
COMPRESSION_SCRIPT = Path(r"E:\20241004_MCI\20260730_Molecular_Neurobiology\05_Data_and_Code\code\compress_mci_centered_panel_20260722.R")
DEVELOPMENT_SCRIPT = Path(r"E:\20241004_MCI\20260730_Molecular_Neurobiology\05_Data_and_Code\code\optimize_mci_development_panel_20260722_161919.py")
DEVELOPMENT_DIR = HISTORICAL_ROOT / "mci_development_optimized_panel_20260722_161919"
DEVELOPMENT_CONSISTENCY = DEVELOPMENT_DIR / "development_gene_consistency.csv"
ARCHIVED_LOCKED_PANEL = DEVELOPMENT_DIR / "development_locked_12_gene_panel.csv"
ARCHIVED_MAPPING_AUDIT = DEVELOPMENT_DIR / "development_locked_12_mapping_audit.csv"
ARCHIVED_QA_NOTES = DEVELOPMENT_DIR / "optimization_QA_notes.txt"
CURRENT_LEDGER = RESULTS / "panel_provenance_ledger.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_source_files() -> list[Path]:
    paths = [
        CANDIDATE_REGISTRY,
        BOOTSTRAP_STABILITY,
        PROVISIONAL_PANELS,
        COMPRESSION_SCRIPT,
        DEVELOPMENT_SCRIPT,
        DEVELOPMENT_CONSISTENCY,
        ARCHIVED_LOCKED_PANEL,
        ARCHIVED_MAPPING_AUDIT,
        ARCHIVED_QA_NOTES,
        CURRENT_LEDGER,
    ]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("required archived selection source is missing: " + "; ".join(missing))
    return paths


def reconstruct_locked_panel(consistency: pd.DataFrame) -> pd.DataFrame:
    required = {
        "gene",
        "trajectory_class",
        "development_consistent",
        "development_min_abs_effect",
    }
    missing = required - set(consistency.columns)
    if missing:
        raise ValueError(f"development consistency table missing required columns: {sorted(missing)}")
    consistent = consistency.loc[consistency["development_consistent"].astype(bool)].copy()
    peak = consistent.loc[consistent["trajectory_class"].eq("MCI_peak")].sort_values(
        "development_min_abs_effect", ascending=False, kind="mergesort"
    ).head(6)
    trough = consistent.loc[consistent["trajectory_class"].eq("MCI_trough")].sort_values(
        "development_min_abs_effect", ascending=False, kind="mergesort"
    ).head(6)
    selected = pd.concat([peak, trough], ignore_index=True)
    selected["selection_rank_within_trajectory"] = selected.groupby("trajectory_class").cumcount() + 1
    selected["selection_data_scope"] = "GSE63060_and_GSE63061_development_only"
    return selected[[
        "gene",
        "trajectory_class",
        "development_min_abs_effect",
        "selection_rank_within_trajectory",
        "selection_data_scope",
    ]]


def audit() -> tuple[pd.DataFrame, dict[str, object]]:
    sources = require_source_files()
    registry = pd.read_csv(CANDIDATE_REGISTRY)
    consistency = pd.read_csv(DEVELOPMENT_CONSISTENCY)
    archived = pd.read_csv(ARCHIVED_LOCKED_PANEL)
    ledger = pd.read_csv(CURRENT_LEDGER)
    selected = reconstruct_locked_panel(consistency)

    archived_symbols = archived["gene"].astype(str).tolist()
    reconstructed_symbols = selected["gene"].astype(str).tolist()
    current_symbols = ledger["source_symbol"].astype(str).tolist()
    selected_counts = selected.groupby("trajectory_class").size().to_dict()
    current_consistent = int(consistency["development_consistent"].astype(bool).sum())

    audit_record: dict[str, object] = {
        "status": "ARCHIVED_LOCKED12_DEVELOPMENT_SELECTION_REPLAY_VERIFIED",
        "candidate_universe_count": int(len(registry)),
        "development_consistent_candidate_count": current_consistent,
        "selected_counts_by_trajectory": {key: int(value) for key, value in selected_counts.items()},
        "selection_rule": (
            "from strict MCI-centered candidate registry, require direction agreement in GSE63060 and "
            "GSE63061; rank each trajectory by minimum absolute age/sex-adjusted effect across both "
            "development cohorts; retain six MCI-peak and six MCI-trough genes"
        ),
        "selection_scope": "GSE63060_and_GSE63061_development_only",
        "gse249477_used_for_selection": False,
        "reconstructed_order_matches_archived_locked_panel": reconstructed_symbols == archived_symbols,
        "reconstructed_symbols_match_current_provenance_ledger": reconstructed_symbols == current_symbols,
        "claim_boundary": (
            "development-only archived replay; not independent discovery, not external validation, "
            "not clinical biomarker validation, and not a reconstruction of the original raw preprocessing lineage"
        ),
        "source_hashes": {str(path): sha256_file(path) for path in sources},
    }
    return selected, audit_record


def main(output_dir: Path = RESULTS) -> int:
    selected, audit_record = audit()
    if not audit_record["reconstructed_order_matches_archived_locked_panel"]:
        raise ValueError("reconstructed archived selection does not match the archived locked panel")
    if not audit_record["reconstructed_symbols_match_current_provenance_ledger"]:
        raise ValueError("reconstructed archived selection does not match the current provenance ledger")
    output_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output_dir / "archived_locked12_selection_audit.csv", index=False)
    (output_dir / "archived_locked12_selection_audit.json").write_text(
        json.dumps(audit_record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit_record, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output_dir))
