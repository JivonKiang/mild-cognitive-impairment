"""Audit inherited cis-eQTL-to-AD MR output against V5 evidence gates."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
LEGACY_RUN = WORKSPACE / "MCI" / "runs" / "Biomni_lab_downloads_20260810_012752" / "mci_paper"
MR_SCRIPT = LEGACY_RUN / "scripts" / "21_mendelian_randomization.R"
COLOC_SCRIPT = LEGACY_RUN / "scripts" / "21b_colocalization.R"
MR_RESULTS = LEGACY_RUN / "results_tables" / "mr_results_cis_eqtl.csv"
COLOC_RESULTS = LEGACY_RUN / "results_tables" / "mr_colocalization.csv"
OUT = ROOT / "03_RESULTS"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def contains(text: str, fragment: str) -> bool:
    return fragment.casefold() in text.casefold()


def audit() -> tuple[list[dict[str, str]], dict[str, object]]:
    mr_text = MR_SCRIPT.read_text(encoding="utf-8")
    coloc_text = COLOC_SCRIPT.read_text(encoding="utf-8")
    with MR_RESULTS.open(encoding="utf-8", newline="") as handle:
        mr_rows = list(csv.DictReader(handle))
    with COLOC_RESULTS.open(encoding="utf-8", newline="") as handle:
        coloc_rows = list(csv.DictReader(handle))

    fdr_positive = sorted({row["gene"] for row in mr_rows if float(row["fdr"]) < 0.05})
    max_pp_h4 = max(float(row["PP.H4"]) for row in coloc_rows)
    findings = [
        {
            "finding_id": "exposure_effect_scale",
            "code_evidence": "beta_exp := Zscore; se_exp := 1",
            "gate_failed": "Recorded exposure beta and standard error",
            "impact": "The inherited output cannot document the actual eQTL effect and uncertainty used for this audit.",
        },
        {
            "finding_id": "instrument_pruning",
            "code_evidence": "distance-based clumping; 1 Mb window",
            "gate_failed": "LD-reference-based instrument clumping",
            "impact": "Distance-only pruning does not satisfy the V5 LD-clumping evidence gate.",
        },
        {
            "finding_id": "coloc_gene_locus_identity",
            "code_evidence": "eqtl_dedup <- eqtl_all[, .SD[which.min(Pvalue)], by = SNP]",
            "gate_failed": "Gene-locus-specific colocalization input",
            "impact": "Global SNP deduplication precedes gene-specific colocalization and does not preserve a gene-locus input record.",
        },
        {
            "finding_id": "outcome_compatibility",
            "code_evidence": "script header: cis-eQTL -> AD; AD GWAS outcome",
            "gate_failed": "MCI-stage phenotype compatibility",
            "impact": "An AD-GWAS outcome cannot establish a causal MCI-stage result in this manuscript.",
        },
        {
            "finding_id": "reproduction_provenance",
            "code_evidence": "legacy script expects /workspace/mci_paper/mr_data inputs; V5 retains output tables but not those raw input files, execution environment or execution log",
            "gate_failed": "Inspectable rerun provenance",
            "impact": "Inherited output tables alone do not meet the V5 reproducibility threshold.",
        },
    ]
    required = {
        "mr_uses_z_as_beta": contains(mr_text, "beta_exp := Zscore"),
        "mr_sets_se_to_one": contains(mr_text, "se_exp := 1"),
        "mr_uses_distance_clumping": contains(mr_text, "distance-based clumping"),
        "mr_targets_ad_outcome": contains(mr_text, "cis-eQTL -> AD"),
        "coloc_globally_deduplicates_eqtls": contains(
            coloc_text, "eqtl_dedup <- eqtl_all[, .SD[which.min(Pvalue)], by = SNP]"
        ),
    }
    if not all(required.values()):
        raise RuntimeError(f"Legacy MR code fingerprint changed: {required}")
    audit_payload: dict[str, object] = {
        "status": "EXCLUDED_FROM_CAUSAL_INTERPRETATION_LEGACY_METHOD_AUDIT",
        "analysis_role": "supplementary code-and-provenance exclusion audit",
        "online_method_source_status": "PENDING_ONLINE_METHOD_SOURCE_VERIFICATION_NETWORK_UNAVAILABLE_2026-08-14",
        "legacy_output_summary": {
            "mr_rows": len(mr_rows),
            "mr_genes": len({row["gene"] for row in mr_rows}),
            "fdr_lt_005_genes": fdr_positive,
            "colocalization_rows": len(coloc_rows),
            "maximum_pp_h4": max_pp_h4,
            "pp_h4_ge_05_count": sum(float(row["PP.H4"]) >= 0.5 for row in coloc_rows),
        },
        "code_fingerprints": required,
        "source_files": {
            str(path): sha256(path)
            for path in [MR_SCRIPT, COLOC_SCRIPT, MR_RESULTS, COLOC_RESULTS]
        },
        "claim_boundary": "The inherited cis-eQTL-to-AD MR and colocalization outputs are excluded from causal, MCI-stage, biomarker and target-validation interpretation. The audit does not judge correctly implemented MR or colocalization methods.",
    }
    return findings, audit_payload


def main() -> int:
    findings, payload = audit()
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "legacy_mr_method_findings.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(findings[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(findings)
    (OUT / "legacy_mr_method_audit.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
