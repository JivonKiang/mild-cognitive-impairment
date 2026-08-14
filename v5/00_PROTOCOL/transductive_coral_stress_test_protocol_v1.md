# Target-Unlabeled CORAL Stress-Test Protocol

Status: exploratory method extension, frozen 2026-08-14 after the original V5
primary analysis.

## Question

Can a batch-level, target-unlabeled covariance-alignment method restore MCI-HC
transport from the AddNeuroMed development family to two public target domains?

This is an exploratory stress test, not a prospective validation design. Its
results cannot alter the pre-existing primary transcriptome result or establish
clinical utility.

## Method source

CORAL is based on second-order source-target alignment without target labels:
Sun, Feng and Saenko, *Return of Frustratingly Easy Domain Adaptation*, AAAI
2016, DOI `10.1609/aaai.v30i1.10306`. The source was verified online on
2026-08-14.

## Cohort roles and analysis unit

- Source: GSE63060 and GSE63061 AddNeuroMed HC/MCI rows only, n=208.
- Target 1: all 587 GSE140829 expression rows before diagnosis metadata access;
  final HC/MCI assessment is 383 sample rows with unresolved subject mapping.
- Target 2: all 62 GSE249477 RNA-seq expression rows before diagnosis parsing;
  final HC/MCI assessment is 41 sample rows.
- Shared feature universe: 11,575 exact symbols present in both source series
  and both target assays. No aliases or target-label-driven feature changes are
  allowed.

## Locked transform and classifier

1. Fit median imputation and `StandardScaler` on the 208 source rows only.
2. Fit randomized PCA on the source-transformed shared-gene matrix only, with
   50 components and seed `20260814`.
3. Transform one complete unlabeled target batch into the source PCA space.
4. CORAL-align source latent features to the target latent mean and covariance
   with ridge regularization `0.001`.
5. Fit a target-unlabeled `StandardScaler` on target latent features; apply it
   to both the aligned source latent matrix and target latent matrix. This
   controls numerical scale after covariance alignment while preserving the
   target-unlabeled batch-level protocol.
6. Fit `LogisticRegression(C=0.2, solver='lbfgs', max_iter=1500)` using aligned
   source latent features and source HC/MCI labels only.
7. Lock all target probabilities. Only then read target diagnosis metadata and
   calculate descriptive AUC, bootstrap 95% CI, Brier score, ECE, and descriptive
   evaluation-set calibration summaries.

## Mandatory exclusion and interpretation rules

- Target diagnosis labels must not be used for imputation, scaling, PCA, CORAL,
  feature intersection, model fitting, hyperparameter selection, thresholding,
  or calibration.
- The procedure uses a full target batch and is therefore a transductive,
  batch-level stress test. It is not a single-sample deployment procedure.
- ComBat is intentionally not run. Although its empirical-Bayes batch-adjustment
  source is verified (Johnson et al., DOI `10.1093/biostatistics/kxj037`), a
  diagnosis-preserving covariate cannot be supplied for the deliberately
  unlabelled target rows. Applying ComBat without that covariate could remove or
  alter the primary MCI-HC signal; reporting it as label-free adaptation would be
  invalid.
- A positive result would still not establish clinical validation. A negative
  result supports only the statement that this specific target-unlabeled
  covariance-alignment approach did not restore transport in these audits.

## Reproducibility and outputs

Run:

```powershell
E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_transductive_coral_audit.py
```

Outputs are `03_RESULTS/transductive_coral_metrics.csv`,
`transductive_coral_predictions_anonymized.csv`, and
`transductive_coral_audit.json`. The audit records the source hashes, target
row counts before label access, method settings, label-access contract, allowed
claims, and forbidden claims.
