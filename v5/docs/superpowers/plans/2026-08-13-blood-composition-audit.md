# Blood Composition Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reproducible two-method blood-cell composition audit without exposing target labels to score generation or modifying the locked MCI diagnostic models.

**Architecture:** A method registry fixes external implementation provenance, a label-blind score generator produces cohort-wide score tables, and a separate post-lock assessor produces FDR-controlled exploratory association and descriptive prediction-score analyses. Existing figure and manuscript infrastructure consume only aggregated, mirrored outputs.

**Tech Stack:** Python 3 bundled runtime, pandas/numpy/scipy/statsmodels where already available, externally installed R packages in an isolated library, unittest, matplotlib.

---

### Task 1: Lock provenance and label-isolation contracts

**Files:**
- Create: `00_PROTOCOL/blood_composition_audit_protocol_v1.md`
- Create: `00_PROTOCOL/blood_composition_method_registry.csv`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing contract tests**

```python
def test_blood_composition_protocol_blocks_target_labels_from_scoring(self):
    text = (ROOT / "00_PROTOCOL" / "blood_composition_audit_protocol_v1.md").read_text(encoding="utf-8")
    self.assertIn("all cohort rows before diagnosis labels", text)
    self.assertIn("never modify, refit, calibrate, or threshold", text)

def test_blood_composition_registry_requires_reproducibility_fields(self):
    rows = list(csv.DictReader((ROOT / "00_PROTOCOL" / "blood_composition_method_registry.csv").open(encoding="utf-8")))
    self.assertEqual({row["method"] for row in rows}, {"xCell", "MCPcounter"})
    for row in rows:
        self.assertTrue(row["source_url"])
        self.assertTrue(row["resolved_commit"])
        self.assertTrue(row["license_status"])
```

- [ ] **Step 2: Run the new tests and verify red**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: failures because the protocol and registry do not yet exist.

- [ ] **Step 3: Add the protocol and registry**

State the exact two-method boundary, source URLs, commit identifiers, accessed
date, license status, executable location and `BLOCKED` fallback rule. Do not
vendor external code or signature files.

- [ ] **Step 4: Re-run the contract tests**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: the two new contract tests pass.

### Task 2: Build a pure label-blind scoring interface

**Files:**
- Create: `02_ANALYSIS/run_blood_composition_audit.py`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing unit tests**

```python
def test_composition_score_generation_is_invariant_to_label_permutation(self):
    module = load_module("blood_composition", ANALYSIS / "run_blood_composition_audit.py")
    expression = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"], index=["r1", "r2"])
    first = module.score_all_rows_label_blind(expression, scorer=lambda x: x.sum(axis=1).to_frame("score"))
    labels = np.array([1, 0])
    second = module.score_all_rows_label_blind(expression, scorer=lambda x: x.sum(axis=1).to_frame("score"))
    self.assertEqual(labels.tolist(), [1, 0])
    pd.testing.assert_frame_equal(first, second)

def test_composition_assessment_rejects_missing_preassessment_rows(self):
    module = load_module("blood_composition_rows", ANALYSIS / "run_blood_composition_audit.py")
    with self.assertRaises(ValueError):
        module.assert_complete_preassessment_rows(["r1", "r2"], ["r1"])
```

- [ ] **Step 2: Run the new tests and verify red**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: failure because the module/functions do not exist.

- [ ] **Step 3: Implement only the tested interface and manifest writer**

Implement `score_all_rows_label_blind`, `assert_complete_preassessment_rows`,
SHA256 input fingerprints and a manifest whose target-label usage fields are
all false before assessment. The script must stop with `BLOCKED` when either
frozen external implementation cannot run; no historical xCell values may be
read.

- [ ] **Step 4: Re-run the contract tests**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: new unit tests pass.

### Task 3: Install frozen external dependencies and execute score generation

**Files:**
- Create: `02_ANALYSIS/install_blood_composition_dependencies.R`
- Create: `03_RESULTS/blood_composition_install_audit.json`
- Create: `03_RESULTS/blood_composition_preassessment_manifest.json`

- [ ] **Step 1: Write a failing audit test for explicit install outcome**

```python
def test_composition_install_audit_has_explicit_method_outcome(self):
    audit = json.loads((ROOT / "03_RESULTS" / "blood_composition_install_audit.json").read_text(encoding="utf-8"))
    self.assertEqual(set(audit["methods"]), {"xCell", "MCPcounter"})
    self.assertTrue(all(row["status"] in {"READY", "BLOCKED"} for row in audit["methods"].values()))
```

- [ ] **Step 2: Run and verify red**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: failure because the install audit has not been created.

- [ ] **Step 3: Install only into an isolated project-local R library and run**

Pin commit sources, capture installer output and session information. Run score
generation on every cohort row before endpoint labels. Emit READY or BLOCKED
for each method. Do not install into global R libraries.

- [ ] **Step 4: Re-run audit test and inspect the manifest**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: status is explicit and no target labels were available during score
generation.

### Task 4: Post-lock exploratory assessment, traceability and figure integration

**Files:**
- Modify: `02_ANALYSIS/run_blood_composition_audit.py`
- Modify: `02_ANALYSIS/make_figure_pack.py`
- Modify: `05_MANUSCRIPT/manuscript_draft_v1.md`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `04_FIGURES/panel_traceability.csv`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests for evidence boundaries**

```python
def test_composition_outputs_describe_scores_not_absolute_proportions(self):
    audit = json.loads((ROOT / "03_RESULTS" / "blood_composition_audit.json").read_text(encoding="utf-8"))
    self.assertTrue(audit["exploratory_only"])
    self.assertFalse(audit["absolute_proportions_reported"])
    self.assertFalse(audit["diagnostic_model_modified"])
```

- [ ] **Step 2: Run and verify red**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py`

Expected: failure because the post-lock audit output does not exist.

- [ ] **Step 3: Implement post-lock tests and integration**

For each score, calculate predeclared MCI-HC effect and within-cohort/method
Benjamini-Hochberg FDR after score lock. Generate a direction concordance table
and prediction-score Spearman correlations with FDR. Add a supplementary
figure only when a method is READY; otherwise add a transparent blocked-status
supplementary table. Mirror all source data and state that findings cannot
rescue transport, establish mechanism or certify a biomarker.

- [ ] **Step 4: Run full QA**

Run: `E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py`

Expected: PASS with zero failures.
