# Post-Lock Composition Source Lock Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the post-lock MCPcounter exploratory analysis reproducible and fail closed when any frozen score, diagnosis-label source, or locked prediction input no longer matches its declared provenance.

**Architecture:** Keep anonymous MCPcounter scores immutable. Before any diagnosis label is accessed, validate the preassessment manifest against the score table and validate the versioned registry against every diagnosis-label source. Verify the pre-existing prediction and audit artifacts independently before descriptive correlations, and report unpooled directional counts explicitly.

**Tech Stack:** Python 3, pandas, scipy, statsmodels, unittest, SHA256.

---

### Task 1: Freeze post-lock input provenance

**Files:**
- Create: `00_PROTOCOL/blood_composition_postlock_source_registry_v1.json`
- Modify: `02_ANALYSIS/run_blood_composition_postlock_assessment.py`
- Test: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing source-lock tests**

Add tests proving that a malformed preassessment manifest fails before labels can be read and that a mismatched official diagnosis-label SHA fails before labels can be read.

- [ ] **Step 2: Run the focused tests to verify failure**

Run: `& E:\AI_Research_Assistant\scripts\utils\python.cmd -m unittest 99_QA.test_contract.ContractTests.test_postlock_manifest_coverage_mismatch_fails_before_label_access 99_QA.test_contract.ContractTests.test_postlock_label_hash_mismatch_fails_before_label_access -v`

Expected: FAIL because the preassessment structure validator and source registry validator do not yet exist.

- [ ] **Step 3: Implement minimal validators**

Validate the manifest row count, column list, cohort set, cohort row counts, score-column counts, and contiguous assessment indices. Load and SHA-validate all declared label sources before the first cohort label loader runs.

- [ ] **Step 4: Run focused tests to verify passing**

Run the Task 1 command again.

### Task 2: Freeze prediction inputs and direction reporting

**Files:**
- Modify: `02_ANALYSIS/run_blood_composition_postlock_assessment.py`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing prediction and direction-summary tests**

Add one test for a prediction SHA mismatch and one test requiring available-effect, MCI-higher, HC-higher, null-direction, and explicit concordance fields.

- [ ] **Step 2: Run focused tests to verify failure**

Run: `& E:\AI_Research_Assistant\scripts\utils\python.cmd -m unittest 99_QA.test_contract.ContractTests.test_postlock_prediction_hash_mismatch_blocks_correlation 99_QA.test_contract.ContractTests.test_postlock_direction_summary_reports_available_and_direction_counts -v`

Expected: FAIL because the locked prediction loader and explicit summary helper do not yet exist.

- [ ] **Step 3: Implement minimal reporting changes**

SHA-validate the frozen predictions and their audit before correlations. Merge only unpooled direction-count summaries into the composition table and record the same summaries in the audit JSON.

- [ ] **Step 4: Run focused tests to verify passing**

Run the Task 2 command again.

### Task 3: Regenerate and release-gate the exploratory outputs

**Files:**
- Create: `03_RESULTS/blood_composition_postlock_summary.csv`
- Create: `03_RESULTS/blood_composition_postlock_prediction_correlations.csv`
- Create: `03_RESULTS/blood_composition_postlock_audit.json`
- Modify: `99_QA/test_contract.py` only if a reproducibility contract is missing

- [ ] **Step 1: Run the post-lock analysis in a disposable output directory**

Run: `& E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_blood_composition_postlock_assessment.py --output-dir <temporary-directory>`

Expected: 40 summary rows and correlations only where exact endpoint mappings are available.

- [ ] **Step 2: Inspect artifacts and run the full QA gate**

Run: `& E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py`

Expected: `PASS` with zero failures.

- [ ] **Step 3: Publish only after all validations pass**

Run the analysis with its canonical `03_RESULTS` output directory. Keep all language and future Figure S1 material explicitly `single-method exploratory`; do not modify primary model results or claims.
