# Cross-Study Stable Signal Audit Design

## Objective

Determine whether a peripheral whole-blood MCI-versus-HC signal that is
direction-concordant in both AddNeuroMed development studies remains
discriminative when the resulting model is locked and transported to the
previously audited external cohorts. This is a reliability analysis, not a
clinical biomarker validation or causal-discovery exercise.

## Scope and units

- Primary estimand: MCI versus HC in peripheral whole-blood bulk mRNA.
- Development family: GSE63060 and GSE63061 only; they are not independent
  external validation cohorts.
- External transport assessments: GSE140829 sample rows and GSE249477 sample
  rows. Their outcome labels are prohibited before locked prediction output.
- All source matrices, labels and platform mapping are already linked by
  `run_pairwise_external.py`; no raw sample identifiers are copied into the
  release package.

## Design

1. Intersect exact-case canonicalized genes across GSE63060 and GSE63061.
2. Generate synchronized, stratified five-fold outer splits within each study.
   At each fold, the held-out partitions from both studies are excluded from
   every selection and model-fit operation.
3. In the paired training partitions, calculate study-specific standardized
   mean differences (Hedges-like pooled-standard-deviation effect estimates)
   for MCI minus HC. Selection requires the same non-zero sign in both studies
   and a prespecified stability score: `min(abs(effect_60), abs(effect_61))`.
   Select the top 100 eligible genes, ordered deterministically by stability,
   then mean absolute effect and symbol. If fewer than 100 are eligible, retain
   all eligible genes and record this condition.
4. Fit a median-imputed, standardized ridge logistic-regression model with
   fixed `C=0.2` on the combined paired training data only. Score held-out rows
   separately for each study. The outer-fold model is compared against a pooled
   nested-feature ridge comparator whose feature selector is fit only on the
   same paired training rows.
5. Lock the stable-gene selection once on all 208 eligible AddNeuroMed
   MCI/HC rows and fit the final ridge model. Score all rows of GSE140829 and
   GSE249477 before their disease labels are read or endpoint rows filtered.
   Target labels are then used only for aggregate MCI-versus-HC AUC, Brier,
   descriptive calibration and bootstrap confidence intervals.
6. Run a development-label permutation null using the identical split,
   selection and modeling path. This tests whether apparent cross-study
   development performance exceeds label-free chance under the full pipeline;
   it does not make a causal claim.

## Required boundaries

- Never call the within-family result independent validation or leave-study-out
  validation.
- Never fit, calibrate, threshold, select, normalize across samples, or adapt
  using external disease labels.
- The external findings stay negative if their confidence intervals include
  chance or their calibration is poor. No method is promoted because it is
  computationally sophisticated.
- The existing raw ridge, rank/z, TSP, conformal, OOD, leakage attack and
  domain-separability audits are retained; this extension is an additional
  stress test, not a replacement.

## Outputs and integration

- A script emits fold-level selection records, OOF predictions, locked
  transport predictions, metrics, and an audit JSON with source hashes and
  explicit target-label isolation booleans.
- Fig. 5 gains a compact panel comparing pooled-ridge and cross-study-stable
  ridge across the paired hold-outs and external assessments. The existing
  TSP/rank analyses remain in the figure or supplement according to available
  space; no existing virtual analyses are removed.
- The methods protocol, manuscript, evidence ledger, code-data-figure map,
  figure manifest/source mirror, README and QA contract are extended together.

## Success criteria

- Tests demonstrate that a discordant gene is excluded, a concordant gene is
  selected deterministically, held-out study rows cannot affect selection, and
  target-label permutation cannot alter predictions before assessment.
- A full run produces no claim boundary violation, and all registered QA and
  figure-source checks pass.
