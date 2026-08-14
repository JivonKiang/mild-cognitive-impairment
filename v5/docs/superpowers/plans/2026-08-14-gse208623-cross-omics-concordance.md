# GSE208623 Cross-Omics Concordance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a source-locked, sex-adjusted public DNA-methylation association and transcriptome-concordance supplement without modifying or validating the existing diagnostic models.

**Architecture:** A resumable official-GEO downloader acquires only the processed beta matrix and official EPIC annotation after Range and size verification. A separate analysis module validates the 60-column sample map, performs MCI-HC CpG association on M values adjusted for sex, then creates aggregate gene-level methylation summaries and descriptive concordance with the locked AddNeuroMed transcriptome effect reference. Every stage emits hashes and blocks rather than substitutes a result when an acquisition or mapping gate fails.

**Tech Stack:** Python 3, urllib.request, gzip, pandas, numpy, scipy, statsmodels, unittest.

---

### Task 1: Freeze protocol and acquisition contracts

**Files:**
- Create: `00_PROTOCOL/analysis_protocol_gse208623_cross_omics_v1.md`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests** requiring the protocol's MCI-HC-only, sex-adjusted, nonpredictive, noncausal and no-IDAT boundaries.
- [ ] **Step 2: Run the selected test** and confirm that it fails before the test-backed module exists.
- [ ] **Step 3: Write the frozen protocol** with expected official file sizes, valid gzip and exact sample-column gates, CpG model, gene mapping, concordance and pathway-block rules.
- [ ] **Step 4: Run the selected test** and confirm it passes.

### Task 2: Implement resumable official acquisition

**Files:**
- Create: `02_ANALYSIS/download_gse208623_official_sources.py`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests** for Content-Range parsing, rejecting a non-206 Range response, rejecting an unexpected complete byte size, and accepting a complete gzip file only after its header can be read.
- [ ] **Step 2: Run the selected tests** and confirm expected failure because the module does not exist.
- [ ] **Step 3: Implement a range-resumable downloader** that stores `.part` files, verifies the official expected byte count, atomically promotes completed files, and writes no analysis output.
- [ ] **Step 4: Run the selected tests** and confirm they pass.

### Task 3: Implement source-gated association and concordance

**Files:**
- Create: `02_ANALYSIS/run_gse208623_cross_omics_audit.py`
- Modify: `99_QA/test_contract.py`

- [ ] **Step 1: Write failing tests** for M-value conversion, sex-adjusted coefficient direction, BH correction, 60-column metadata alignment, rejection of invalid beta values, and no classifier or causal-analysis API.
- [ ] **Step 2: Run the selected tests** and confirm expected failure.
- [ ] **Step 3: Implement the module** with source hashes, MCI-HC-only OLS models, unambiguous official gene mapping, deterministic aggregate gene Z summaries and descriptive transcriptome-concordance summaries.
- [ ] **Step 4: Run selected tests** and confirm they pass.

### Task 4: Publish and integrate only after input gates pass

**Files:**
- Create: `03_RESULTS/gse208623_cross_omics_cpg_association.csv`
- Create: `03_RESULTS/gse208623_cross_omics_gene_concordance.csv`
- Create: `03_RESULTS/gse208623_cross_omics_audit.json`
- Modify: `99_QA/run_full_qa.py`
- Modify: `05_MANUSCRIPT/claim_evidence_ledger.md`
- Modify: `05_MANUSCRIPT/supplementary_analysis_index.md`
- Modify: `README.md`

- [ ] **Step 1: Write failing publication and regeneration tests** for source hashes, aggregate-only outputs, no model-fitting flags and byte-identical regeneration.
- [ ] **Step 2: Run the selected tests** and confirm expected failure before a complete official source acquisition.
- [ ] **Step 3: Download, validate and analyze the sources; if any source gate fails, publish only an explicit blocked audit and leave primary claims unchanged.**
- [ ] **Step 4: Run the contract suite and full QA** and inspect fresh exit status before reporting a result.
