# Top-Scoring Pairs Transport Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` task-by-task. This package is not a Git repository; the approved project-guard check is the write gate.

**Goal:** Test whether a locked, development-only relative-expression Top-Scoring Pairs (TSP) classifier transports across GSE140829 microarray and GSE249477 RNA-seq public data without target-domain fitting or target-label use before assessment.

**Architecture:** Within each AddNeuroMed development fold, rank genes by `f_classif`, retain the fixed top 100 genes, select the ten highest-separation disjoint gene pairs using only development labels, and score each sample by the proportion of MCI-oriented pair orderings. Fit no target-domain normalizer, classifier, calibration model, threshold, or feature-selection operation. Produce OOF development results and one locked final assessment for each external cohort; target labels are used only to calculate endpoint metrics after predictions are fixed.

**Tech Stack:** Python 3.10, pandas, NumPy, scikit-learn, official GEO matrices, unittest.

---

### Task 1: Lock TSP invariants with failing tests

**Files:**
- Modify: `99_QA/test_contract.py`
- Create: `02_ANALYSIS/run_tsp_relative_expression_transport.py`

- [ ] **Step 1: Add failing behavior tests**

```python
def test_tsp_pair_selection_is_development_label_driven_and_disjoint(self):
    tsp = load_module("tsp_relative_expression_transport", ANALYSIS / "run_tsp_relative_expression_transport.py")
    x = pd.DataFrame(
        {"A": [5, 4, 1, 2], "B": [1, 2, 5, 4], "C": [4, 5, 2, 1], "D": [2, 1, 4, 5]}
    )
    pairs = tsp.select_tsp_pairs(x, np.array([1, 1, 0, 0]), candidate_gene_count=4, pair_count=2)
    self.assertEqual(len(pairs), 2)
    self.assertEqual(len({gene for pair in pairs for gene in pair[:2]}), 4)

def test_tsp_target_predictions_do_not_depend_on_target_labels(self):
    tsp = load_module("tsp_relative_expression_transport_target", ANALYSIS / "run_tsp_relative_expression_transport.py")
    pairs = [("A", "B", 1), ("C", "D", -1)]
    target = pd.DataFrame({"A": [4, 1], "B": [1, 4], "C": [4, 1], "D": [1, 4]})
    np.testing.assert_allclose(tsp.predict_pair_vote(target, pairs), tsp.predict_pair_vote(target, pairs))
```

- [ ] **Step 2: Verify the tests fail**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
```

Expected: failure because `run_tsp_relative_expression_transport.py` does not exist.

- [ ] **Step 3: Implement only the tested pair selection and vote routines**

`select_tsp_pairs(expression, y, candidate_gene_count=100, pair_count=10)` must:

1. use `f_classif` only on the supplied development expression and labels;
2. retain at most 100 finite-score genes with deterministic lexical tie-breaking;
3. compute the absolute MCI-versus-HC ordering-probability gap for every candidate pair;
4. orient each retained pair so a positive vote is the ordering more frequent in MCI;
5. greedily retain up to ten pairs without reusing a gene; and
6. return deterministic `(gene_high, gene_low, orientation)` tuples where `orientation=1` means `gene_high > gene_low` votes MCI and `-1` means the reverse ordering votes MCI.

`predict_pair_vote(expression, pairs)` must return the mean oriented vote across pairs and must accept only expression and frozen pairs.

- [ ] **Step 4: Verify the new tests pass**

Run the same command. Expected: tests pass.

### Task 2: Generate locked internal and two-cohort transport assessments

**Files:**
- Modify: `02_ANALYSIS/run_tsp_relative_expression_transport.py`
- Generate: `03_RESULTS/tsp_transport_metrics.csv`
- Generate: `03_RESULTS/tsp_transport_selected_pairs.csv`
- Generate: `03_RESULTS/tsp_transport_predictions_anonymized.csv`
- Generate: `03_RESULTS/tsp_transport_audit.json`

- [ ] **Step 1: Add a failing artifact and target-label-boundary test**

```python
def test_tsp_transport_audit_is_label_free_before_external_assessment(self):
    audit = json.loads((ROOT / "03_RESULTS" / "tsp_transport_audit.json").read_text(encoding="utf-8"))
    self.assertEqual(audit["status"], "RELATIVE_EXPRESSION_TRANSPORT_AUDIT")
    self.assertFalse(audit["gse140829_target_labels_used_for_pair_selection"])
    self.assertFalse(audit["gse249477_target_labels_used_for_pair_selection"])
    self.assertTrue(audit["gse140829_target_labels_used_only_for_assessment"])
    self.assertTrue(audit["gse249477_target_labels_used_only_for_assessment"])
    self.assertEqual(audit["fixed_pair_count"], 10)
```

- [ ] **Step 2: Verify it fails, then implement**

Use `load_dev` and `load_external` from `run_pairwise_external.py` and
`load_target_expression` / `parse_assessment_labels` from
`run_gse249477_rank_transport.py`.

For the primary endpoint only:

1. create AddNeuroMed HC/MCI development data from the exact shared-gene intersection;
2. obtain five-fold OOF predictions by selecting pairs independently in every training fold;
3. select the final ten pairs only from all 208 AddNeuroMed HC/MCI rows;
4. run pair-vote inference over every GSE140829 retained expression row before using its HC/MCI labels for metrics;
5. rank/vote every GSE249477 target row before calling `parse_assessment_labels`, then retain HC/MCI only for its locked metrics;
6. use `metric_row` and `bootstrap_auc` for aggregate discrimination/calibration descriptions, but never fit target calibration or thresholds; and
7. export only aggregate metrics, gene pairs, sequential `assessment_index` predictions, hashes, and flags.

The metric rows must be exactly:

```text
NESTED_INTERNAL_SENSITIVITY,GSE63060_GSE63061
PRIMARY_EXTERNAL_TRANSPORT_AUDIT,GSE140829
SECONDARY_INDEPENDENT_CROSS_PLATFORM_TRANSPORT_AUDIT,GSE249477
```

- [ ] **Step 3: Regenerate and check output**

Run:

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_tsp_relative_expression_transport.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
```

Expected: artifact contract passes; GSE249477 is still assessment-only and all 62 target rows are predicted before label parsing.

### Task 3: Make regeneration and source data mandatory

**Files:**
- Modify: `99_QA/run_full_qa.py`
- Modify: `02_ANALYSIS/make_figure_pack.py`

- [ ] **Step 1: Add failing full-QA checks**

Require the four TSP result artifacts, byte-identical regeneration in a temporary `--output-dir`, source-script SHA-256 binding, three expected statuses, two external assessment counts (383 and 41), ten final disjoint pairs, no raw target identifiers, and false target-label pair-selection flags.

- [ ] **Step 2: Register the script and source-data mirrors**

Register `run_tsp_relative_expression_transport.py` in both script lists.
Mirror metrics, selected pairs and audit JSON to `04_FIGURES/source_data`; retain prediction rows only in `03_RESULTS` unless a later figure requires them.

- [ ] **Step 3: Run full QA**

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py
```

Expected: zero failures and byte-identical TSP regeneration.

### Task 4: Integrate the evidence without promoting it to validation

**Files:**
- Modify: `02_ANALYSIS/make_figure_pack.py`
- Modify: `05_MANUSCRIPT/manuscript_draft_v1.md`
- Modify: `05_MANUSCRIPT/figure_captions_v3.md`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `05_MANUSCRIPT/source_verification_register.md`
- Modify: `05_MANUSCRIPT/Code_Data_Figure_Mapping.md`
- Modify: `05_MANUSCRIPT/submission_readiness_v1.md`
- Modify: `README.md`

- [ ] **Step 1: Add TSP to Fig. 5 as a method comparison**

Fig. 5a compares raw transcriptome, rank/z, and TSP AUCs for GSE140829 and GSE249477 where defined. The title and caption must say `development-only relative-expression TSP`, `transport audit`, and `not external validation`.

- [ ] **Step 2: Record literature and scope**

Record official PubMed URLs and access date for `PMID 16646797` and
`PMID 23496976`. State that the records support a relative-expression/TSP
method family for cross-platform molecular classification; they do not support
MCI efficacy, clinical utility, or any expectation of positive performance.

- [ ] **Step 3: Rebuild and verify**

```text
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\make_figure_pack.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\test_contract.py
E:\AI_Research_Assistant\scripts\utils\python.cmd 99_QA\run_full_qa.py
```

Expected: six nonblank figures, complete source-data hashes, and zero QA failures.
