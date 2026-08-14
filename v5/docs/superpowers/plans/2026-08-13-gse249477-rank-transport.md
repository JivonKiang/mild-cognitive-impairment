# GSE249477 Rank Transport Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` task-by-task. This package is not a Git repository, so the project-guard check substitutes for a branch/worktree gate.

**Goal:** Add a reproducible, label-leakage-free transport audit of the existing rank/z-score transcriptome sensitivity in the independent GSE249477 whole-blood RNA-seq cohort.

**Architecture:** The audit derives per-sample ranks using all 62 expression rows before any endpoint filtering, intersects exact gene symbols with the AddNeuroMed development data, and fits feature selection and logistic regression only within AddNeuroMed folds. GSE249477 diagnosis labels are accessed after all target-domain transformations and only to compute a single locked assessment. The result is a secondary independent cross-platform boundary audit, never clinical validation or a replacement primary result.

**Tech Stack:** Python 3.10, pandas, NumPy, scikit-learn, official GEO series-matrix title metadata, official GSE249477 TPM matrix, unittest.

---

### Task 1: Lock the target-label and assay boundaries

**Files:**
- Modify: `99_QA/test_contract.py`
- Create: `02_ANALYSIS/run_gse249477_rank_transport.py`

- [ ] **Step 1: Write the failing contract**

```python
def test_gse249477_rank_transport_is_target_label_free_before_assessment(self):
    audit_path = ROOT / "03_RESULTS" / "gse249477_rank_transport_audit.json"
    self.assertTrue(audit_path.exists())
    result = json.loads(audit_path.read_text(encoding="utf-8"))
    self.assertEqual(result["status"], "SECONDARY_INDEPENDENT_TRANSPORT_AUDIT")
    self.assertEqual(result["assessment_cohort"], "GSE249477")
    self.assertEqual(result["assessment_unit"], "sample-row; 21 HC and 20 MCI")
    self.assertFalse(result["target_labels_used_for_transform"])
    self.assertFalse(result["target_labels_used_for_fit"])
    self.assertFalse(result["target_labels_used_for_feature_selection"])
    self.assertFalse(result["target_labels_used_for_threshold_or_calibration"])
    self.assertTrue(result["target_labels_used_only_for_assessment"])
```

- [ ] **Step 2: Verify the contract fails**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
```

Expected: `FAIL` because `gse249477_rank_transport_audit.json` does not yet exist.

- [ ] **Step 3: Implement the smallest self-contained transport script**

Create `02_ANALYSIS/run_gse249477_rank_transport.py` with these concrete rules:

```python
GSE249477_RAW = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_raw_count_normalize_04-10-2025.csv.gz"
)
GSE249477_SERIES = Path(
    r"E:\AI_Research_Assistant\data\operations\run_bundles\mci_gse249477_task9\GSE249477_series_matrix.txt.gz"
)

def sample_rank_transform(frame: pd.DataFrame) -> pd.DataFrame:
    ranks = np.apply_along_axis(
        lambda row: pd.Series(row).rank(method="average").to_numpy(), 1, frame.to_numpy(dtype=float)
    )
    return pd.DataFrame(
        (ranks - np.nanmean(ranks, axis=1, keepdims=True))
        / np.where(np.nanstd(ranks, axis=1, keepdims=True) > 0, np.nanstd(ranks, axis=1, keepdims=True), 1.0),
        index=frame.index,
        columns=frame.columns,
    ).astype(np.float32)
```

The script must:

1. Parse all 62 public series-matrix titles into `HC`, `MCI`, or `AD` after confirming one unique title-derived identifier per TPM column.
2. Canonicalize symbols by case only, collapse duplicate RNA-seq gene rows by median, and use the exact intersection with the two AddNeuroMed development matrices.
3. Transform all 62 GSE249477 rows sample-wise before restricting to the 21 HC and 20 MCI assessment rows.
4. Use five-fold AddNeuroMed-only nested evaluation, then one final AddNeuroMed fit. Imputation, `SelectKBest(f_classif, k=min(200, feature_count))`, scaling and ridge logistic regression with fixed `C=0.2` must remain inside the development-only pipeline.
5. Write only aggregate metrics, anonymized assessment indices, source SHA-256 hashes and explicit allowed/forbidden claims to `03_RESULTS`; never export raw target expression or sample identifiers.

- [ ] **Step 4: Generate the audit and rerun contracts**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_gse249477_rank_transport.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
```

Expected: all contracts pass and the audit states that GSE249477 labels were used only for assessment.

### Task 2: Enforce reproducibility and figure-source mirroring

**Files:**
- Modify: `99_QA/run_full_qa.py`
- Modify: `02_ANALYSIS/make_figure_pack.py`
- Generate: `03_RESULTS/gse249477_rank_transport_metrics.csv`
- Generate: `03_RESULTS/gse249477_rank_transport_predictions_anonymized.csv`
- Generate: `03_RESULTS/gse249477_rank_transport_audit.json`
- Generate: `04_FIGURES/source_data/gse249477_rank_transport_metrics.csv`
- Generate: `04_FIGURES/source_data/gse249477_rank_transport_audit.json`

- [ ] **Step 1: Add a failing full-QA assertion**

Add checks that require all three results artifacts, require `SECONDARY_INDEPENDENT_TRANSPORT_AUDIT`, exactly 41 `HC`/`MCI` assessment rows, and reject any audit with a true target-label transform, fit, feature-selection, calibration or threshold flag.

```python
check(len(gse249_rank_metrics) == 2, "GSE249477 rank transport must contain nested and external rows", failures)
check(set(gse249_rank_metrics["status"]) == {"NESTED_INTERNAL_SENSITIVITY", "SECONDARY_INDEPENDENT_TRANSPORT_AUDIT"}, "GSE249477 rank transport statuses changed", failures)
check(int(gse249_rank_audit["assessment_label_counts"]["HC"]) == 21 and int(gse249_rank_audit["assessment_label_counts"]["MCI"]) == 20, "GSE249477 rank transport label count changed", failures)
```

- [ ] **Step 2: Verify full QA fails before regeneration**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py
```

Expected: `FAIL` because the new script/results are not yet registered.

- [ ] **Step 3: Register the script and source-data mirrors**

Add the transport script to `model_scripts` and `scripts` in `run_full_qa.py`. Add the metrics CSV and audit JSON to `copy_source_data()` in `make_figure_pack.py`. Do not add per-sample prediction CSV to public figure source data because it has no panel consumer; retain it only under `03_RESULTS` for reproducibility.

- [ ] **Step 4: Rebuild figures and run full QA**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\make_figure_pack.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py
```

Expected: six figures remain nonblank, manifest hashes refresh, and full QA reports zero failures.

### Task 3: Record complete candidate discovery and preserve claim boundaries

**Files:**
- Create: `01_DATA_AUDIT/public_mci_blood_candidate_discovery_20260813.md`
- Modify: `05_MANUSCRIPT/source_verification_register.md`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `05_MANUSCRIPT/manuscript_draft_v1.md`
- Modify: `05_MANUSCRIPT/submission_readiness_v1.md`
- Modify: `README.md`

- [ ] **Step 1: Record the exact public discovery set**

The discovery note must name both retained NCBI GEO raw files:

```text
01_DATA_AUDIT/gds_mild_cognitive_impairment_blood_search_20260813.json
01_DATA_AUDIT/gds_mild_cognitive_impairment_blood_summary_20260813.json
```

It must report 38 exact `mild cognitive impairment` plus blood GEO series metadata records on 2026-08-13, and state that no previously unregistered independent whole-blood bulk mRNA MCI-HC cohort was eligible. It must identify GSE249477, GSE270454, GSE140829 and AddNeuroMed as the only compatible records in this discovery set, with their distinct roles.

- [ ] **Step 2: Add precise manuscript wording**

Allowed wording:

```text
A prespecified sample-wise rank/z-score sensitivity was independently assessed in GSE249477 after target-domain transformation without diagnosis labels; its labels were used only to quantify the locked transport result.
```

Forbidden wording:

```text
GSE249477 clinically validates the model; the rank transform solves domain shift; or results prove a generalizable MCI blood biomarker.
```

The manuscript must report the observed AUC and bootstrap interval only after the script has generated them. It must state the 41-row sample count and retain the cross-platform RNA-seq designation.

- [ ] **Step 3: Verify final package**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py
```

Expected: contracts and full QA pass; source-data hashes in `04_FIGURES/figure_manifest.json` match regenerated files.
