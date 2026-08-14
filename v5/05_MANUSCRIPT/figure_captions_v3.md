# Figure Captions and Statistical Contract (V5)

All AUC intervals are bootstrap 95% intervals over the stated analysis unit.
For GSE140829, the unit is explicitly a sample row because a validated
sample-to-subject map is unavailable. Logistic calibration intercept and slope
are descriptive assessments fit on evaluation labels after predictions were
locked; they are not external recalibration parameters. ECE uses ten fixed
probability bins unless otherwise stated. Supplementary/stress-test results are
not positive validation claims.

## Fig. 1 | Cohort admission and evidence boundaries

The candidate-compatibility registry admits the AddNeuroMed development family
to the primary peripheral whole-blood bulk-mRNA MCI-versus-HC estimand, followed
by a GSE140829 external sample-row transport audit and GSE249477 secondary
cross-platform transport stress test. GSE282742 progression and GSE208623
cross-omics analyses are separate, non-pooled branches. PBMC, intervention,
cell-specific, brain and non-coding-RNA sources are excluded from primary model
performance. Counts are sample rows unless a unique-subject count is explicitly
shown. The figure is a design and evidence-boundary schematic, not a validation,
mechanistic or causal result.

## Fig. 2 | Leakage and transportability

(a) Repeated five-fold diagnostic comparison of fold-contained selection and a
full-development label-aware selection attack. (b) Bidirectional transfer
between the two AddNeuroMed series and transport to GSE140829. (c) Primary
MCI-versus-HC transcriptome AUC in development and GSE140829. Points show AUC
and bars show bootstrap 95% intervals. The external result is a sample-row
transport audit with high OOD rate, not a positive clinical validation or
subject-level validation.

## Fig. 3 | Primary-endpoint external uncertainty

All panels are restricted to GSE140829 `stage_mci_vs_hc` rows (n=383; HC=249,
MCI=134). (a) Five display-bin descriptive calibration curves for the age-plus-
sex and transcriptome models. ECE is computed separately with ten fixed
probability bins, as stated in the metrics table. Per-bin sample counts are
annotated; no binomial uncertainty interval is implied by the connected points.
(b) Training-derived standardized distance by label with the development
95th-percentile threshold; this distance is a diagnostic quantity, not a
validated domain-membership probability. (c) Selective error against external
coverage using thresholds indexed only by development confidence quantiles;
external labels are assessment-only.
The prediction rows and selective-risk table are included in figure source
data.

## Fig. 4 | Exact-symbol mapping boundary

The frozen signed 12-gene score maps exactly to 12/12 development symbols but
only 10/12 symbols in GSE140829 (LAT1-3TM and ZSWIM6 absent), GSE249477 (PHF15
and LAT1-3TM absent), and GSE270454 (PHF15 and LAT1-3TM absent). For GSE270454,
the count-column to GEO metadata join was independently audited as unique for
45/45 whole-blood RNA samples before age/sex adjustment. External partial scores
are sensitivity or boundary analyses; missing symbols are not substituted. AUC
intervals are bootstrap 95% intervals over sample rows.

## Fig. 5 | Development-family signal exceeds a label null but does not restore external transport

(a) A label-free sample-wise rank/z-score transform gives AUC 0.605 in the
GSE140829 sample-row audit and AUC 0.576 (bootstrap 95% CI 0.387-0.776) in the
secondary independent cross-platform GSE249477 transport audit (n=41; HC=21,
MCI=20). For GSE249477, all 62 target rows were transformed before diagnosis
labels were parsed; labels were used only for the locked final HC/MCI assessment.
No target calibration was fit. The GSE249477 Brier score was 0.488 and ECE was
0.488, so it is a weak and uncertain ranking stress test with severe probability
miscalibration. The same forest plot reports a development-only
relative-expression top-scoring-pairs (TSP) transport audit. Ten disjoint,
MCI-oriented pairs were selected only from AddNeuroMed. All 587 GSE140829 rows
and all 62 GSE249477 rows were scored before endpoint assessment; TSP AUCs were
0.548 (95% CI 0.486-0.612; n=383) and 0.465 (95% CI 0.300-0.629; n=41),
respectively. The same locked forest includes target-unlabeled latent CORAL: all
587 GSE140829 and 62 GSE249477 target expression rows were available as batches
before diagnosis labels were accessed, but labels were not used for alignment,
model selection or prediction. AUCs remained 0.544 (95% CI 0.484-0.597; n=383)
and 0.421 (95% CI 0.255-0.610; n=41), respectively. CORAL is a batch-level
transductive stress test, not a single-sample deployment model or external
validation. (b) GSE140829 AUC, ECE and Brier for raw, rank/z and CORAL scores
are shown in metric-specific axes and should not be interpreted on a common
clinical utility scale. (c)
Development-calibrated split-conformal prediction at nominal coverage 0.90 is
evaluated only in the GSE140829 external sample-row audit; observed mean
coverage is 0.669 across 20 repeats. The whisker shows the 2.5th-97.5th
percentile across repeats. Rank/z and TSP are method stress tests, not external
validation, clinical deployment, or evidence that relative-expression ordering
solves domain shift. (d) The observed synchronized AddNeuroMed development
out-of-fold AUC for the cross-study stable ridge is compared with 100
cohort-preserving label permutations, each repeating feature selection and
outer-fold evaluation. The empirical P value tests development-family
discrimination only; it is not a test of external validation, clinical utility
or deployment validity.

The lower panels add a synchronized cross-study stable-signal ridge model and
its label-null reference distribution.
Its mapping-only 11,575-gene candidate universe was fixed before target disease
labels were accessed. In GSE63060 and GSE63061 outer folds, direction-concordant
MCI-HC effects were selected only from paired training partitions; the resulting
stable model gave AUC 0.840 and 0.779, respectively. A 100-repeat development
label-permutation path gave p=0.0099 for the pooled development OOF AUC. After
final selection on all AddNeuroMed MCI/HC rows, all 587 GSE140829 and all 62
GSE249477 rows were scored before endpoint assessment. Locked external AUCs
were 0.541 (95% CI 0.479-0.598; n=383) and 0.524 (95% CI 0.500-0.575; n=41).
This establishes development-family directional consistency but not independent
validation, calibrated transport, or a solution to domain shift.

## Fig. 6 | Phenotype-metadata-independent domain separability during model fitting

The domain classifier distinguishes AddNeuroMed development sample rows from
GSE140829 external sample rows using expression values and a binary domain
label. All rows are retained for the primary domain audit; disease labels are
not used for preprocessing, feature selection, thresholding or domain-model
fitting. They are used for row alignment and post hoc phenotype-stratified
sensitivities, which are reported in Table S6.
(a) Receiver-operating-characteristic curve; (b) predicted probability
distributions by domain; (c) AUC across ten repeated 5-fold cross-validation
runs. The aggregate domain AUC is 1.000 (bootstrap 95% CI 1.000-1.000;
permutation p=0.00020). This is evidence of domain separation, not a disease
classifier, and GSE140829 remains a sample-row audit because subject
independence is unresolved.

## Fig. 7 | Cross-cohort MCI-HC effect concordance and descriptive heterogeneity

The fixed 11,575-gene exact-symbol universe was constructed before GSE140829
or GSE249477 diagnosis labels were accessed. (a) Pairwise Spearman
correlations of independently estimated MCI-minus-HC Hedges' g. (b) Fraction
of nonzero effects with the same direction; GSE249477 comparisons use 10,923
estimable genes because 652 shared genes are non-estimable in that cohort.
(c) Cochran Q/I-squared distribution restricted to GSE63060 and GSE63061,
reported descriptively without a heterogeneity P value or released pooled
effect. GSE140829 is a sample-row assessment with unresolved subject
independence and is excluded from pooled inference. The figure quantifies
cohort-level concordance or discordance; it does not establish biological
mechanism, cell-composition causality, clinical utility or validation.

## Fig. S1 | Source-locked GSE208623 cross-omics association audit

GSE208623 is an independent-source Chinese peripheral-blood-leukocyte EPIC
cohort. The sex-adjusted MCI-minus-HC association uses 20 MCI and 20 HC sample
rows; 20 AD rows are excluded before modeling. (a) Accounting of all 738,242
processed beta-matrix CpGs, the CpGs with Benjamini-Hochberg FDR below 0.05,
and genes with an unambiguous official EPIC symbol mapping. (b) Aggregate
signed methylation Z statistics versus the frozen mean GSE63060/GSE63061
MCI-minus-HC transcriptome Hedges' g in the 10,422 exact-symbol overlap; the
annotated Spearman rho is descriptive. (c) Same-sign and opposite-sign
fractions among sign-eligible genes, with the mean from 10,000 fixed
gene-sign permutations and the corresponding empirical P value.

This is a source-locked, sex-adjusted cross-omics association audit, not
prediction, mRNA-model validation, replication, methylation-expression
mediation, cell-composition proof, biological mechanism, Mendelian
randomization, or causal inference. Age and blood-cell composition are not
adjusted. Gene-level signed-Z values aggregate CpG associations and are not
independent gene-level differential-methylation tests.
