# Cross-Cohort Effect Heterogeneity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a source-locked, non-predictive audit of MCI-HC gene-effect concordance and heterogeneity across the four compatible public blood-transcriptome cohorts.

**Architecture:** The module reuses existing expression loaders to build a mapping-only shared exact-symbol universe without diagnosis metadata for target construction. It then computes each cohort's independent standardized effects, descriptive pairwise concordance, and two-study AddNeuroMed heterogeneity while blocking all cross-cohort pooled inference.

**Tech Stack:** Python, pandas, NumPy, SciPy, existing project loaders, unittest.

---

### Task 1: Add failing contracts and lock the protocol

**Files:**
- Create: `00_PROTOCOL/analysis_protocol_cross_cohort_effect_heterogeneity_v1.md`
- Modify: `99_QA/test_contract.py`

- [ ] Add tests for Hedges' g sign, finite sampling variance, pairwise concordance, and rejection of unmatched gene universes.
- [ ] Run selected tests and confirm missing-module failure.
- [ ] Confirm protocol forbids pooled four-cohort inference and target-label-driven selection.

### Task 2: Implement the minimal effect audit

**Files:**
- Create: `02_ANALYSIS/run_cross_cohort_effect_heterogeneity.py`
- Modify: `99_QA/test_contract.py`

- [ ] Compute independent Hedges' g and variance per cohort in a predeclared universe.
- [ ] Compute pairwise Spearman/sign summaries and development-only two-study Q/I-squared.
- [ ] Emit aggregate CSV/JSON artifacts with source hashes and label-use audit fields.
- [ ] Run the selected contracts to green.

### Task 3: Publish and integrate release gates

**Files:**
- Create: `03_RESULTS/cross_cohort_effects.csv`
- Create: `03_RESULTS/cross_cohort_effect_pairwise_summary.csv`
- Create: `03_RESULTS/cross_cohort_effect_heterogeneity.csv`
- Create: `03_RESULTS/cross_cohort_effect_audit.json`
- Modify: `99_QA/run_full_qa.py`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `05_MANUSCRIPT/supplementary_analysis_index.md`

- [ ] Add regeneration and byte-equality checks.
- [ ] Run the script to publish artifacts.
- [ ] Run contract tests and full QA before reporting conclusions.
