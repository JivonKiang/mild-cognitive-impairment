# MCI Public Reliability Pairwise Protocol v1

Status: FROZEN BEFORE RESULT GENERATION.

This protocol is a new, isolated protocol. The prior v2 audit package remains
unchanged and is treated as historical evidence. No result from this package
may be described as a clinical diagnostic validation without an independent
sample-level validation design.

## 1. Endpoints

- Primary endpoint: `stage_mci_vs_hc`.
- Key secondary endpoint: `stage_mci_vs_ad`.
- Context endpoint: `stage_ad_vs_hc_context`.
- Progression endpoint: `progression_p_mci_vs_s_mci`; it is not a stage
  endpoint and cannot validate the primary endpoint.

The primary analysis uses one development study family (GSE63060/GSE63061,
AddNeuroMed) and one independent external assessment (GSE140829). Because the
two AddNeuroMed series share a study family, they are not counted as independent
external domains. GSE140829 is an external transport assessment; it is not
called a positive validation result unless subject independence and all other
compatibility requirements are established.

## 2. Public data and source boundaries

Development inputs:

- GSE63060 official matrix and official sample-label table.
- GSE63061 official matrix and official sample-label table.

External input:

- GSE140829 complete GEO series matrix with 587 samples.
- GPL15988 official platform table for probe-to-symbol mapping.
- GSE249477 official count/TPM matrix with 62 samples (HC=21, MCI=20,
  AD=21) for a cross-platform stress test only; it is not a primary external
  validation cohort.

Only public, source-linked data are allowed. Raw matrices and sample-level
identifiers remain outside this output directory. Results store source paths,
source hashes, aggregate counts, and execution hashes.

## 3. Frozen biological comparators

The locked 12-gene stage score is retained exactly as a comparator:

`GPI, PHF15, VAT1, EIF3B, UCP2, LAT1-3TM, C22ORF39, USP38, ZSWIM6, PRDM10, RP2, GCA`.

Peak genes are weighted +1/6 and trough genes -1/6. Platform symbols may be
canonicalized by case only (for example `C22orf39` to `C22ORF39`); this is not
an alias substitution. No aliases, gene substitutions, 12+7 union panel, or
post-result re-selection is allowed.

If a target platform lacks an exact symbol, the full locked score is BLOCKED
for that cohort. A partial 10-gene score may be reported only as an explicitly
labeled sensitivity/audit result and may never be promoted to the primary
model. In GSE140829, LAT1-3TM and ZSWIM6 are currently absent from the official
GPL15988 symbol mapping, so the full external locked score remains blocked.

## 4. Leakage-free transcriptome model

The new model is a separate, prespecified computational comparator intended
to use the public data efficiently without information leakage.

- Input genes are the exact-symbol intersection after source-platform mapping.
- Probe aggregation is unsupervised and performed without target labels.
- For every development fold, feature screening, imputation, scaling, and
  model fitting occur inside the training fold. The ridge penalty `C=0.2` is
  fixed before execution and is not selected after inspecting results.
- Feature screening uses `SelectKBest(f_classif, k=min(200, n_features))` inside
  the fold, followed by ridge logistic regression.
- The external GSE140829 labels are never used for fitting, feature selection,
  scaling, calibration, threshold choice, OOD threshold choice, or adaptation.
- The final external model is fit once on all eligible development samples
  after internal nested evaluation is complete.

The age-only, sex-only, and age-plus-sex models are prespecified baselines and
are all executed in the primary and study-held-out analyses.
The transcriptome model is not automatically promoted if it fails to add
value over age-plus-sex or if its external result is unstable.

## 5. Validation and uncertainty

Development stability is evaluated with stratified folds and a bidirectional
study-held-out check between GSE63060 and GSE63061. Independent external
assessment is reported separately for GSE140829. GSE249477 is reported only as
a cross-platform stress test. Required metrics are AUC, Brier score,
descriptive logistic-calibration intercept and slope, ECE with bin count, and
subject/sample bootstrap intervals. Calibration intercept and slope are
evaluation-set descriptive assessments fit after predictions are locked; they
do not modify external predictions and must not be described as external
recalibration. The method is recorded explicitly; it is not an OLS regression
of the binary outcome on the logit score.

Selective prediction is reported as a transportability analysis only:

- OOD score is based on training-derived standardized feature distance.
- The abstention threshold is fixed from development-only predictions.
- Risk-coverage is evaluated on the external cohort using thresholds indexed by
  development-only confidence quantiles; external labels are never used to
  choose the threshold.

An additional OOD sensitivity audit compares the primary mean absolute
standardized distance with training-only Ledoit-Wolf Mahalanobis and robust PCA
distances. Thresholds are fixed from development q90, q95 and q99 values.
External labels are used only for post hoc error-rate and discrimination
assessment; they are not used to fit a distance model or choose a threshold.

Prespecified sensitivity analyses additionally test (i) a sample-wise rank and
z-score transform that uses no cohort-level target labels and (ii) split
conformal prediction with a development-only calibration split. External
labels are used only to assess discrimination, calibration, coverage and
singleton error; they are never used to tune the transformation, conformal
quantile, threshold or model.

The cross-study stable-signal audit is a distinct sensitivity analysis. Its
candidate universe is restricted before target diagnosis metadata are accessed
to exact symbols present in both AddNeuroMed development studies, GSE140829 and
GSE249477. In synchronized five-fold outer partitions of GSE63060 and GSE63061,
MCI-minus-HC standardized effects are computed separately in the two training
partitions. Only non-zero direction-concordant effects are eligible; the top
100 are ordered by the minimum absolute study effect, then mean absolute effect
and symbol. A fixed C=0.2 ridge model is fit on the combined paired training
rows and scored separately on both held-out studies. A pooled training-only
`SelectKBest` ridge comparator uses the same mapping universe. Final feature
selection occurs once on all 208 AddNeuroMed MCI/HC rows, after which every
GSE140829 and GSE249477 row is scored before target disease labels are read for
assessment. A 100-repeat development-label permutation null re-runs stable
feature selection and the outer-fold model; it evaluates only development-family
discrimination. Neither development-family result is independent validation,
and targets are never used for selection, preprocessing, fitting, prediction
calibration, thresholding or adaptation. After predictions are locked, target
labels are used only for descriptive performance and calibration assessment.

No clinical utility, treatment benefit, causal effect, or prospective
performance claim is permitted.

## 6. Required interpretation

The primary manuscript story is: label leakage can inflate small-cohort
transcriptomic claims, and a frozen, leakage-aware public-data workflow can
separate apparent discovery performance, same-family transfer, independent
transportability and uncertainty. The MCI case study is an application of that
reliability framework, not a clinical deployment study.

MR, virtual knockout, LINCS, network analyses, docking, MD, and ADME remain
supplementary hypothesis-generating evidence. They cannot rescue a blocked
expression mapping or create external validation.
