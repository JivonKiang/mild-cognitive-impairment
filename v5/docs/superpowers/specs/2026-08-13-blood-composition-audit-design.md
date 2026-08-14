# Blood Composition Audit Design

## Purpose

Add a public-data-only, leakage-aware blood-cell composition audit to explain
whether transcriptome-model scores may covary with blood-cell abundance
signals. It is supplementary evidence for the MCI-versus-HC reliability and
transportability manuscript, not a diagnostic-model modification or biological
validation experiment.

## Preconditions and evidence boundary

- Primary estimand remains peripheral whole-blood bulk mRNA, MCI versus HC.
- GSE63060 and GSE63061 are one AddNeuroMed development family; GSE140829 is a
  sample-row external transport assessment; GSE249477 is a secondary RNA-seq
  cross-platform audit.
- Every cohort row is scored before diagnosis labels are read, filtered, or
  supplied to any cell-score procedure.
- Cell scores are method-specific transcriptomic enrichment/abundance scores,
  never absolute cell percentages or ground truth.
- The locked diagnostic models, feature sets, score calibration and thresholds
  are immutable. Model-score versus cell-score analysis is descriptive only.
- Diagnosis association is exploratory, performed only after score lock,
  corrected within cohort across predeclared cell categories, and never used to
  select a score method or alter any model.

## Method options considered

### A. xCell only

Uses the published xCell enrichment score method (Aran et al., Genome Biology,
2017; PMID 29141660). It is familiar and broad, but the repository snapshot
does not expose a top-level license and it has more dependencies. It is not
enough for an independent technical check.

### B. xCell plus MCP-counter-style coarse-score concordance (selected)

Runs two independent published score families on exactly the same untouched
cohort rows. xCell provides broad enrichment scores; MCP-counter provides
coarser lineage-oriented abundance scores (Becht et al., Genome Biology, 2016;
PMID 27765066). Only predeclared coarse categories with a documented mapping
between methods are tested. A conclusion requires directional consistency, not
just significance in one method. This is the best available public-data-only
option if both frozen implementations can run reproducibly.

### Controlled MCPcounter recovery

The current R environment cannot fetch CRAN packages or GitHub archives because
of a TLS credential failure. The MCPcounter fixed repository commit remains
accessible through the GitHub API, including the GPL-3 package DESCRIPTION,
the unmodified source function and the CC0 signature table. If package
installation remains blocked, the audit may source those exact, hash-verified
upstream files only in a disposable temporary directory and invoke the upstream
`appendSignatures` function. This is not a substitute signature method. xCell
remains `BLOCKED`, and all MCPcounter results are therefore single-method
exploratory evidence rather than the selected two-method concordance result.

### C. CIBERSORTx-style deconvolution

Would estimate fractions but entails a separate web service/reference-matrix
dependency and less auditable versioning in this package. It is excluded from
the current implementation.

## Architecture and data flow

1. A version-and-license registry captures URLs, resolved commit hashes,
   licenses, accessed date, frozen inputs and non-redistribution boundary for
   each external implementation.
2. A label-blind scoring wrapper loads each whole cohort expression matrix,
   maps each platform to exact HGNC symbols, records row/gene counts and emits
   scores before calling an endpoint-label parser.
3. A post-lock assessment module reads labels solely for predeclared MCI-HC
   tests, Benjamini-Hochberg FDR within cohort/method and effect-direction
   meta-summary. It uses no test result for model or method selection.
4. A descriptive model-score module joins existing anonymized locked
   predictions to the already-locked composition scores by row order and
   computes rank correlations with multiplicity correction. It never changes a
   model object or its predictions.
5. Source data, audit JSON, supplementary table and a composition supplementary
   figure are mirrored into the figure traceability system. Language explicitly
   limits results to confounding/context evidence.

## Failure handling

- A missing implementation, absent licence clarity, failed installation, or
  incompatible platform mapping leaves that method `BLOCKED`; it is not
  replaced by historical result files or unverified gene sets. The only
  permitted recovery is a hash-verified invocation of MCPcounter's unmodified
  upstream source function and CC0 signature file in a disposable temporary
  directory; it does not unblock xCell or create a second method.
- A cohort failing score generation is reported with its cause and omitted from
  cross-cohort summary rather than silently dropped.
- One-method-only or inconsistent direction results are reported as discordant
  exploratory findings.

## Acceptance criteria

- Tests prove all-row preassessment scoring is invariant to a target-label
  permutation and that endpoint labels cannot reach score generation.
- Tests prove composition analyses do not invoke diagnostic model fitting,
  feature selection, calibration or thresholding.
- Tests require method/version/URL/license audit fields and explicit score-scale
  wording.
- Full package QA passes after integration. Any method blocked by reproducibility
  constraints remains explicitly blocked in the manuscript and readiness
  ledger.
