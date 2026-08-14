# Blood Composition Post-Lock Assessment v1

Status: PRESPECIFIED AFTER LABEL-BLIND SCORE LOCK AND BEFORE DIAGNOSIS ASSOCIATION.

## Preconditions and scope

This module may read diagnosis labels only after verifying the SHA256 and all-row
coverage of `03_RESULTS/blood_composition_preassessment_scores_anonymized.csv`
against its companion preassessment manifest. It must then verify every
diagnosis-label source and any locked-prediction artifact against
`blood_composition_postlock_source_registry_v1.json`, before the corresponding
source is read. It consumes the locked MCPcounter scores and must not regenerate,
transform, select, calibrate, threshold, refit, or otherwise alter diagnostic
models or composition scores.

The available evidence is a single-method exploratory MCPcounter source-bridge
audit. xCell remains BLOCKED. Scores are method-specific whole-blood
transcriptomic scores, not absolute cell proportions, flow-cytometry truth, or
immune-mechanism measurements.

## Analysis unit and endpoint

The analysis unit is one public-data sample row. GSE140829 is explicitly a
sample-row audit with unresolved subject independence. Within each cohort,
the endpoint is MCI versus HC; AD and unlabelled rows are excluded only after
the locked all-row score table is verified. No cohorts are pooled for a common
effect estimate or inferential test.

## Predeclared score association analysis

For every MCPcounter score in each cohort, report MCI-minus-HC mean difference,
Hedges' g with the small-sample correction, and a two-sided Welch t-test.
The Benjamini-Hochberg procedure adjusts the ten cell-score p values within
each cohort and MCPcounter method. Report exact unadjusted and adjusted p
values, sample sizes, and effect direction. These are exploratory summaries;
no normality claim is made from the test, and no score is dropped based on its
result.

Cross-cohort output reports only the number of available cohort effects and
the counts with MCI-higher, HC-higher, and no-mean-difference directions, plus
their explicit concordance state. Do not pool cohort effects, calculate a
meta-analytic p value, or make a cross-cohort biological conclusion.

## Locked-prediction correlations

A descriptive Spearman correlation may be computed only for a pre-existing,
locked model prediction table after validating one-to-one agreement of cohort,
assessment index, and binary MCI/HC labels with the post-lock assessment table.
The correlation family is adjusted with Benjamini-Hochberg within the
cohort/model pair. A failed mapping produces `NOT_EVALUATED_MAPPING_MISMATCH`,
not a substituted model or remapped score. These correlations are descriptive,
not mediation, confounding correction, mechanism, or model improvement tests.

## Reporting boundary

The results and any figure must say `single-method exploratory`. They must not
claim cell composition, infiltration, absolute abundance, cross-method
concordance, causal mediation, biomarker validity, external validation,
clinical utility, or recovery of failed external transport. Null and discordant
findings remain in the supplementary record.
