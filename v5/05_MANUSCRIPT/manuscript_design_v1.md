# Manuscript design v5: A reliability audit of public blood-transcriptome MCI models

## Central question

How much of an apparently strong MCI blood-transcriptome result survives strict leakage control, study-held-out transfer, independent external assessment, cross-platform stress testing, and explicit domain-shift diagnostics, after candidate public cohorts are screened by assay compatibility rather than sample count alone?

## Main story

The manuscript is a reliability/transportability study with an MCI case study. It does not claim discovery of a clinically validated diagnostic panel. The primary endpoint is MCI versus HC. MCI versus AD is secondary. P-MCI versus S-MCI is a separate exploratory progression branch. Newly discovered public series are screened in a compatibility registry: GSE18309 is excluded because of its very small PBMC cohort, and GSE312139 is retained only as an orthogonal cell-response context because its readout is not direct blood expression.

## Main figures

1. Design and cohort hierarchy: frozen analysis, two development series, independent microarray external cohort, RNA-seq stress test, and separate progression branch.
2. Leakage and transportability: repeated leakage attack versus clean nested selection; bidirectional study-held-out transfer; independent external collapse and OOD rate.
3. External uncertainty: calibration, training-derived domain distance, and risk-coverage/selective prediction.
4. Locked-panel boundary: exact mapping fractions and partial-score performance across internal, GSE140829 and GSE249477.
5. Domain-invariant sensitivity: label-free sample-wise rank/z-score, a development-only relative-expression TSP transport audit, and a direction-concordant cross-study stable-signal audit test whether invariant feature constructions recover external discrimination; none stably restores transport, while calibration and conformal coverage remain unstable.
6. Unlabelled domain separability: a domain-only classifier distinguishes the AddNeuroMed development domain from GSE140829 before disease labels are used, quantifying the transport boundary.
7. Cross-cohort effect concordance: cohort-wise MCI-minus-HC effects show related-study agreement but cross-platform discordance, with two-study descriptive heterogeneity confined to the AddNeuroMed family.

## Methods commitments

- All target labels are excluded from preprocessing, feature selection, calibration, threshold choice and adaptation for external assessment.
- `SelectKBest` is fit inside training folds; ridge logistic regression uses fixed `C=0.2`.
- The exact 12-gene panel is not replaced, expanded or alias-substituted.
- Full external score is blocked when an exact symbol is absent; partial scores are sensitivity analyses only.
- GSE140829 is the independent external transport assessment; subject-level
  independence remains unresolved and is reported as a sample-row audit.
- GSE249477 is a secondary independent cross-platform transport audit, not
  positive replication or external validation; its target labels are
  assessment-only and its rank/z probabilities are not target calibrated.
- The development-only relative-expression TSP extension selects ten disjoint
  pairs only in AddNeuroMed. It scores all GSE140829 rows and all GSE249477
  rows before endpoint assessment; it fits no target calibration, threshold or
  adaptation and is not external validation or clinical deployment.
- The cross-study stable-signal audit uses a mapping-only candidate universe
  available in both development studies and both target platforms. Inside
  synchronized development-family outer folds it selects only direction-
  concordant training effects, then locks a 100-gene ridge model. Its
  development-family results are not independent validation; GSE140829 and
  GSE249477 diagnosis labels remain assessment-only after all target rows are
  scored.
- Age-only, sex-only and age-plus-sex baselines are reported alongside the
  transcriptome model.
- Candidate public cohorts are screened before analysis using `00_PROTOCOL/candidate_compatibility_registry.csv`; incompatible assays are not pooled into the primary estimand.
- Progression uses subject-level units and is not used as primary-stage validation.
- Domain separability uses a binary domain label (AddNeuroMed development versus GSE140829 external); all rows are retained for the primary domain audit and phenotype metadata are discarded before domain preprocessing, feature selection, thresholding and model fitting. HC/MCI/AD-stratified sensitivities are post hoc only.
- The frozen-panel provenance ledger distinguishes historical source-table direction signs from newly trained coefficients and records exact-symbol boundaries per cohort.
- GSE108215 is a supplementary, paired MCI intervention-response audit with
  35 traceable pre/post subjects and no HC group; it cannot validate the stage
  endpoint or support an intervention-efficacy claim. GSE153881 is blocked
  because official GEO metadata and supplementary filenames conflict on MCI/AD
  group assignments.

## Results order

1. Cohort ledger and endpoint separation.
2. Leakage attack: modest but measurable internal optimism.
3. Study-family transfer versus independent transport.
4. External OOD, calibration and development-threshold selective prediction.
5. Frozen-panel mapping boundary.
6. Label-free rank/z sensitivity, development-only relative-expression TSP
   transport audit, cross-study stable-signal audit, and conformal failure
   boundary.
7. Unlabelled domain separability.
8. Cross-cohort effect concordance and within-family descriptive heterogeneity.
9. Exploratory progression, paired intervention-response context, and
supplementary computational hypotheses.

## Submission positioning

The strongest defensible contribution is methodological: public-data biomarker
claims should report leakage attacks, study-family structure, exact mapping,
external transport, clinical baselines, calibration and OOD. A clinical
biomarker or causal-mechanism framing is not supported by the current evidence.
