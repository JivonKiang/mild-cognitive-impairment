from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "03_RESULTS"
PANEL_SOURCE = Path(r"E:\20241004_MCI\20260730_Molecular_Neurobiology\05_Data_and_Code\source_data\source_figure2_locked_12_gene_panel.csv")
CONTRAST_SOURCE = Path(r"E:\20241004_MCI\20260730_Molecular_Neurobiology\05_Data_and_Code\source_data\source_figure3_locked_score_contrasts.csv")
LEGACY_RESULTS = Path(r"E:\20260805 Biomni\MCI_INTEGRATION_20260810\06_FROZEN_SCORE\frozen_score_results.csv")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mapping_statuses() -> dict[str, dict[str, str]]:
    panel = pd.read_csv(PANEL_SOURCE)
    genes = panel["gene"].astype(str).tolist()
    statuses = {gene: {} for gene in genes}
    for gene in genes:
        statuses[gene]["GSE63060_GSE63061"] = "exact_symbol_present"
    gse140 = pd.read_csv(RESULTS / "metrics.csv")
    missing_140 = set(
        str(gse140.loc[
            (gse140["cohort"] == "GSE140829")
            & (gse140["model"] == "frozen_12_gene_score")
            & (gse140["status"] == "PARTIAL_EXTERNAL_SENSITIVITY"),
            "missing_locked_symbols",
        ].iloc[0]).split(";")
    )
    for gene in genes:
        statuses[gene]["GSE140829"] = "missing_exact_symbol" if gene in missing_140 else "exact_symbol_present"
    stress = json.loads((RESULTS / "gse249477_stress_test.json").read_text(encoding="utf-8"))
    missing_249 = set(stress["missing_locked_symbols"])
    for gene in genes:
        statuses[gene]["GSE249477"] = "missing_exact_symbol" if gene in missing_249 else "exact_symbol_present"
    audit = json.loads((RESULTS / "gse270454_boundary_audit.json").read_text(encoding="utf-8"))
    missing_270 = set(audit["frozen_panel_mapping"]["missing_genes"].split(";"))
    for gene in genes:
        statuses[gene]["GSE270454"] = "missing_exact_symbol" if gene in missing_270 else "exact_symbol_present"
    legacy = pd.read_csv(LEGACY_RESULTS)
    exact = legacy[legacy["analysis"].astype(str) == "exact_mapping_GSE282742"]
    for gene in genes:
        values = exact.loc[exact["metric"].astype(str) == gene, "value"]
        statuses[gene]["GSE282742"] = "exact_symbol_present" if len(values) and values.iloc[0] != "" else "unresolved_exact_symbol"
    return statuses


def main() -> int:
    panel = pd.read_csv(PANEL_SOURCE)
    statuses = mapping_statuses()
    panel_hash = sha256_file(PANEL_SOURCE)
    contrast_hash = sha256_file(CONTRAST_SOURCE)
    rows: list[dict[str, object]] = []
    for _, row in panel.iterrows():
        gene = str(row["gene"])
        trajectory = str(row["trajectory_class"])
        sign = 1 if trajectory == "MCI_peak" else -1
        rows.append(
            {
                "source_symbol": gene,
                "official_symbol": gene,
                "official_symbol_status": "source token retained; independent nomenclature verification not asserted",
                "trajectory_class": trajectory,
                "fixed_direction_sign": sign,
                "fixed_weight": sign / 6.0,
                "source_table": str(PANEL_SOURCE),
                "source_table_sha256": panel_hash,
                "contrast_table": str(CONTRAST_SOURCE),
                "contrast_table_sha256": contrast_hash,
                "lock_status": "frozen legacy panel; no reselection or weight change in V5",
                "selected_before_current_evaluation": True,
                "selection_scope": "archived deterministic GSE63060/GSE63061 development-only selection replay verified; not independent discovery or external validation",
                "exact_alias_refseq_rule": "exact source symbol only; aliases and RefSeq substitution not accepted",
                "platform_probe_mapping_rule": "development exact symbol matrix; GSE140829 official GPL15988 probe annotation; GSE249477 official gene-level matrix",
                "multi_probe_rule": "mean across duplicate exact-symbol probes/rows when present; no target-label probe selection",
                "GSE63060_GSE63061_mapping": statuses[gene]["GSE63060_GSE63061"],
                "GSE140829_mapping": statuses[gene]["GSE140829"],
                "GSE249477_mapping": statuses[gene]["GSE249477"],
                "GSE270454_mapping": statuses[gene]["GSE270454"],
                "GSE282742_mapping": statuses[gene]["GSE282742"],
                "missing_symbol_boundary": ";".join(
                    dataset for dataset, status in statuses[gene].items() if status != "exact_symbol_present"
                ) or "none",
            }
        )
    out = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(RESULTS / "panel_provenance_ledger.csv", index=False)
    audit = {
        "status": "PROVENANCE_LEDGER_GENERATED",
        "panel_size": len(out),
        "source_panel_sha256": panel_hash,
        "source_contrast_sha256": contrast_hash,
        "historical_selection_claim": "archived deterministic development-selection replay is separately verified in archived_locked12_selection_audit.json; this remains development-only, not independent discovery or external validation",
        "alias_substitution_allowed": False,
        "weight_rule": "trajectory direction sign divided by 6; this is a fixed audit comparator, not a newly trained coefficient vector",
        "source_files": [str(PANEL_SOURCE), str(CONTRAST_SOURCE), str(LEGACY_RESULTS)],
        "source_hashes": {str(path): sha256_file(path) for path in [PANEL_SOURCE, CONTRAST_SOURCE, LEGACY_RESULTS]},
    }
    (RESULTS / "panel_provenance_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
