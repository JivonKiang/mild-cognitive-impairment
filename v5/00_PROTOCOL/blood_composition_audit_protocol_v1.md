# Blood Composition Audit Protocol v1

Status: PRESPECIFIED BEFORE COMPOSITION-RESULT GENERATION.

## Purpose and estimand

This public-data-only supplementary audit tests whether whole-blood
transcriptome model scores covary with published transcriptome-derived blood
cell abundance/enrichment scores. The primary estimand remains peripheral
whole-blood bulk mRNA, MCI versus HC. These are not absolute cell percentages;
the audit does not estimate absolute cell percentages, establish an immune
mechanism, or validate a biomarker.

## Cohorts and roles

- GSE63060 and GSE63061: AddNeuroMed development-family cohorts.
- GSE140829: external sample-row transport assessment; subject independence is
  unresolved.
- GSE249477: secondary RNA-seq cross-platform transport audit.

No cohort is added to the primary endpoint solely because composition scores
can be generated.

## Methods and implementation control

The audit first attempts to run fixed snapshots of xCell and MCPcounter in a
project-local R library only. The source URL, resolved commit, package version,
license field, installation output, and R session information are captured in
the method registry and installation audit. External package code and signature
resources are not copied into this package. Historical xCell output files are
prohibited as inputs.

If project-local package installation is unavailable but a fixed MCPcounter
source commit and its CC0 `genes.txt` signature blob can be retrieved from the
official GitHub API, a controlled source-bridge recovery is permitted. It must
write the unmodified upstream `MCPcounter.R` blob and `genes.txt` blob only to
a disposable system temporary directory, verify both API blob SHA values, and
invoke the upstream `appendSignatures` function directly. The bridge must not
reimplement, edit, or persist the MCPcounter function or signatures in this
project. This is still one MCPcounter method, not a new marker method.

This recovery cannot unblock xCell. A failed xCell dependency installation
leaves xCell `BLOCKED`. Accordingly, a source-bridge result supports only a
single-method exploratory audit and cannot support a two-method concordance
conclusion.

Only predeclared coarse lineage-level cell categories with documented
cross-method comparability are summarized. A cell score is a method-specific
transcriptomic enrichment/abundance score, not an absolute proportion or a
cytometry measurement.

## Label isolation and model immutability

For each method and cohort, all cohort rows before diagnosis labels are read,
parsed, filtered, or supplied to the scoring procedure are provided to the
scorer. A complete preassessment row manifest is written before any diagnosis
association is calculated. Target labels are not used for platform mapping, score
generation, preprocessing, method selection, score calibration, thresholding,
or adaptation.

The pre-existing MCI diagnostic models and their predictions are locked. This
audit must never modify, refit, calibrate, or threshold a diagnostic model.
Associations between locked predictions and composition scores are descriptive
rank correlations performed only after both values are locked.

## Post-lock exploratory assessment

Only after all rows have scored, explicit MCI/HC diagnosis labels may be read
to calculate predeclared within-cohort score effects and two-sided tests. The
Benjamini-Hochberg procedure controls FDR within each cohort and method.
Cross-cohort summaries report effect direction and method concordance, not a
pooled biological effect or a validated mechanism.

Any association, correlation, or multi-cohort directional pattern is
exploratory. Null, discordant, unavailable, or one-method-only findings are
reported without selective deletion. The audit cannot rescue poor external
transport, create external validation, prove causality, or support clinical
deployment.

## Failure and reporting rules

If a package fails to install, a source snapshot is unavailable, a required
platform mapping fails, or a scorer errors, the respective method/cohort is
marked `BLOCKED` with a reason. It is not replaced with a different gene set,
historical result, or web-service output. A composition figure is generated
only from `READY` methods; otherwise the supplementary record reports the
blocked state and preserves the primary claims unchanged.
