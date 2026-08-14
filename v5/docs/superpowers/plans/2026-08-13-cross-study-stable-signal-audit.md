# Cross-Study Stable Signal Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a leakage-safe cross-study stable disease-signal audit to the public MCI reliability package.

**Architecture:** A focused analysis script will calculate direction-concordant effects from paired GSE63060/GSE63061 training partitions, train a fixed ridge model, and produce only locked external assessments. Existing figures and manuscript artifacts consume the new aggregate tables without changing the established primary model.

**Tech Stack:** Python, pandas, NumPy, scikit-learn, matplotlib, existing package QA.

---

### Task 1: Specify the behavior with failing contract tests

**Files:**
- Modify: `99_QA/test_contract.py`
- Test: `99_QA/test_contract.py`

- [ ] **Step 1: Add an isolated selector test**

```python
stable = load_module("stable_signal", ANALYSIS / "run_cross_study_stable_signal.py")
effects = pd.DataFrame({"GSE63060": [0.8, 0.9, -0.4], "GSE63061": [0.3, -0.1, -0.5]}, index=["CONCORDANT", "DISCORDANT", "NEGATIVE"])
selected = stable.select_concordant_genes(effects, top_k=10)
self.assertEqual(selected["gene"].tolist(), ["CONCORDANT", "NEGATIVE"])
```

- [ ] **Step 2: Run the focused test and observe the missing-module failure**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py ContractTests.test_cross_study_selector_excludes_direction_discordance`

Expected: failure because `run_cross_study_stable_signal.py` does not yet exist.

- [ ] **Step 3: Add target-label isolation and audited-output tests**

```python
self.assertFalse(audit["gse140829_target_labels_used_for_selection"])
self.assertFalse(audit["gse249477_target_labels_used_for_fit"])
self.assertTrue(audit["gse140829_target_labels_used_only_for_assessment"])
self.assertGreaterEqual(len(metrics), 4)
```

- [ ] **Step 4: Run the focused tests and confirm they fail for absent artifacts**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: only the newly added stable-signal tests fail because the feature is absent.

### Task 2: Implement the leakage-safe stable-signal audit

**Files:**
- Create: `02_ANALYSIS/run_cross_study_stable_signal.py`
- Modify: `00_PROTOCOL/analysis_protocol_pairwise_v1.md`
- Modify: `README.md`

- [ ] **Step 1: Implement pure, testable selection functions**

```python
def study_effects(expression: pd.DataFrame, y: np.ndarray) -> pd.Series:
    """Return MCI-minus-HC standardized effects using only supplied training rows."""

def select_concordant_genes(effects: pd.DataFrame, top_k: int = 100) -> pd.DataFrame:
    """Keep only non-zero same-sign effects and deterministically rank stability."""
```

- [ ] **Step 2: Implement paired outer folds and pooled ridge comparator**

Use `StratifiedKFold(n_splits=5, shuffle=True, random_state=20260813)` separately within each development study. For each synchronized fold, calculate effects from the two training partitions, fit median-imputation plus `StandardScaler` plus `LogisticRegression(C=0.2, penalty="l2", solver="lbfgs", max_iter=1500)`, and store study-specific held-out predictions.

- [ ] **Step 3: Implement final locked external scoring and permutation null**

Load all external expression rows before parsing endpoint labels; map only source/platform information, apply the frozen selected-gene set, write prediction rows, then parse MCI/HC labels solely for metrics. Repeat the paired-fold path under a deterministic development-label permutation and write a null-distribution table.

- [ ] **Step 4: Emit auditable results**

Write `cross_study_stable_metrics.csv`, `cross_study_stable_predictions_anonymized.csv`, `cross_study_stable_selected_genes.csv`, `cross_study_stable_selection_events.csv`, `cross_study_stable_permutation_null.csv`, and `cross_study_stable_audit.json`. Include source hashes, feature counts, comparison identity and every target-label isolation invariant.

- [ ] **Step 5: Run new analysis and contract tests**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_cross_study_stable_signal.py`

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: results exist and all contract tests pass.

### Task 3: Integrate figure, manuscript, and traceability

**Files:**
- Modify: `02_ANALYSIS/make_figure_pack.py`
- Modify: `05_MANUSCRIPT/manuscript_draft_v1.md`
- Modify: `05_MANUSCRIPT/manuscript_design_v1.md`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `05_MANUSCRIPT/Code_Data_Figure_Mapping.md`
- Modify: `05_MANUSCRIPT/submission_readiness_v1.md`

- [ ] **Step 1: Add source-data mirroring and a Figure 5 stable-signal panel**

Render pooled ridge versus cross-study stable ridge with AUC bootstrap intervals for both within-family held-outs and both locked external assessments. Use the exact audit language `development-family cross-study consistency` and `external transport audit`.

- [ ] **Step 2: Update prose without changing evidence tiers**

Describe the analysis as a directional-consistency stress test. State whether it does or does not rescue external transport using its actual results. Retain virtual knockout, MR, LINCS, docking, MD and ADME as supplementary hypothesis generation.

- [ ] **Step 3: Regenerate figures and review traceability**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\make_figure_pack.py`

Expected: six figures, refreshed source mirrors and a current hash manifest.

### Task 4: Full verification

**Files:**
- Modify: `99_QA/run_full_qa.py`
- Test: `99_QA/run_full_qa.py`

- [ ] **Step 1: Register the analysis and assert key artifact invariants**

Check that the final audit has target-label isolation booleans, the selector has no direction-discordant genes, and the figure/source mirror includes each stable-signal result table.

- [ ] **Step 2: Run the complete reproducibility sequence**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_cross_study_stable_signal.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\make_figure_pack.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py
```

Expected: contract suite and full QA report zero failures.

- [ ] **Step 3: Inspect Fig. 5 at publication scale**

Confirm panel labels, confidence intervals, sample-size text and evidence-boundary notes remain readable and do not overlap.
