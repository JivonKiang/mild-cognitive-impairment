# GSE282742 Progression Association Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a deterministic, subject-level, covariate-adjusted exploratory progression-association supplement for GSE282742 without training or evaluating a prediction model.

**Architecture:** The module resolves one baseline row per P-MCI/S-MCI subject, checks exact count and metadata contracts, then applies a log2(counts-per-million + small offset) linear model for each gene with age and sex covariates. A distinct paired P-MCI-to-AD sensitivity uses within-subject log-expression changes. Outputs contain aggregate gene/pathway summaries and a provenance audit, never participant identifiers or predictive scores.

**Tech Stack:** Python 3, pandas, numpy, scipy, statsmodels, gzip, unittest.

---

### Task 1: Lock the protocol and analysis-unit tests

**Files:**
- Create: `00_PROTOCOL/analysis_protocol_gse282742_progression_association_v1.md`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests** for the required protocol wording, exactly 52 unique baseline P-MCI/S-MCI subjects, a single excluded missing-age row, and exactly 11 paired P-MCI-to-AD subjects.
- [ ] **Step 2: Run the selected tests** and verify failure because the analysis module and protocol do not exist.
- [ ] **Step 3: Write the protocol** fixing analysis unit, covariates, transformation, multiplicity correction, and forbidden predictive/causal claims.
- [ ] **Step 4: Run the protocol test** and verify it passes.

### Task 2: Implement deterministic association and paired sensitivity

**Files:**
- Create: `02_ANALYSIS/run_gse282742_progression_association.py`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests** for duplicate-subject rejection, age/sex complete-case handling, paired transition restriction, source hashing, deterministic output, and absence of classifier-fitting APIs.
- [ ] **Step 2: Run the selected tests** and verify expected failure due to the missing module.
- [ ] **Step 3: Implement minimal loaders and statistical functions:** resolve the count orientation and exact sample map, compute log2 CPM, fit ordinary least-squares per gene with P-MCI indicator plus age/sex, calculate Benjamini-Hochberg FDR, and compute paired P-MCI-to-AD changes.
- [ ] **Step 4: Run selected tests** and verify they pass.

### Task 3: Add source-locked pathway and cross-branch summaries

**Files:**
- Create: `00_PROTOCOL/gse282742_progression_gene_set_registry_v1.json`
- Modify: `02_ANALYSIS/run_gse282742_progression_association.py`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests** requiring a pinned gene-set registry, a pathway table derived only from baseline adjusted gene statistics, and explicit exploratory/non-causal labels.
- [ ] **Step 2: Run tests** and verify failure.
- [ ] **Step 3: Implement a small registry-pinned Hallmark/MSigDB-derived pathway summary plus directional overlap with the pre-existing 100-gene stable set; no model fitting is introduced.
- [ ] **Step 4: Run selected tests** and verify they pass.

### Task 4: Publish, regenerate, and integrate QA/manuscript boundaries

**Files:**
- Create: `03_RESULTS/gse282742_progression_baseline_association.csv`
- Create: `03_RESULTS/gse282742_progression_paired_transition.csv`
- Create: `03_RESULTS/gse282742_progression_pathway_summary.csv`
- Create: `03_RESULTS/gse282742_progression_audit.json`
- Modify: `99_QA/run_full_qa.py`
- Modify: `00_PROTOCOL/domain_registry.csv`
- Modify: `01_DATA_AUDIT/data_ledger.md`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `05_MANUSCRIPT/supplementary_analysis_index.md`
- Modify: `05_MANUSCRIPT/manuscript_draft_v1.md`
- Modify: `README.md`

- [ ] **Step 1: Write failing QA assertions** for fresh regeneration and byte-equality of all four outputs and their claim boundaries.
- [ ] **Step 2: Run QA assertions** and verify failure before publication.
- [ ] **Step 3: Generate the four public artifacts, mirror only aggregate source data when appropriate, and update project documentation without changing the primary transport conclusion.**
- [ ] **Step 4: Run contract tests and full QA**; inspect exact exit status and output report before claiming completion.
