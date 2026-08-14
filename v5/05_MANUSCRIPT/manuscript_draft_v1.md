# Leakage-aware transportability audits reveal a reliability boundary for public blood-transcriptome models of mild cognitive impairment

## Manuscript status

This is an evidence-bound first draft for a Nature-family or high-impact medical
journal submission. It is a reliability and transportability study, not a
clinical diagnostic validation report. Quantitative claims are bound to the V5
source files listed in `05_MANUSCRIPT/Code_Data_Figure_Mapping.md`. Background
and citations are classified in `source_verification_register.md` as verified,
pending final matching, or hypothesis-generation sources. No acceptance,
impact-factor or journal-quartile outcome is guaranteed.

## One-sentence argument

In public blood-transcriptome modelling of mild cognitive impairment, we show
that fold-contained feature selection can separate modest internal optimism from
transportability failure, supported by study-held-out transfer, an independent
sample-row audit, domain-shift and uncertainty analyses, with the boundary that
the current evidence does not establish a clinically validated biomarker.

## Abstract

Blood transcriptome models are increasingly proposed for molecular assessment of
mild cognitive impairment (MCI), yet their apparent performance can be affected
by label leakage, study-family structure and shifts in assay or population
composition. We developed a public-data reliability audit that separates these
failure modes rather than pooling incompatible cohorts into a single estimate.
The primary endpoint was MCI versus healthy control (HC). Features were selected,
preprocessed and fitted within development folds; a repeated label-aware
selection attack, bidirectional study-held-out transfer, an independent external
transport assessment, a cross-platform stress test, calibration, out-of-
distribution (OOD) sensitivity and split-conformal prediction were evaluated.

Within the AddNeuroMed development family, the nested transcriptome comparator
achieved an AUC of 0.795 (95% CI 0.722-0.860; n = 208). A full-scope
label-aware selection attack gave an AUC of 0.822 versus 0.802 for the clean
nested diagnostic, corresponding to a mean optimism of 0.0198. Bidirectional
transfer between the two AddNeuroMed series yielded AUCs of 0.770 and 0.730.
In the independent GSE140829 sample-row transport audit, the primary AUC was
0.483 (95% CI 0.420-0.542; n = 383), with a Brier score of 0.344, an expected
calibration error of 0.340 and a development q95 OOD rate of 0.859. A label-free
within-sample rank/z-score transform partially improved external ranking (AUC
0.605; 95% CI 0.538-0.663), but did not establish calibrated transportability.
The separately locked rank/z sensitivity was then assessed as a secondary
independent cross-platform transport audit in GSE249477. After all 62 target
rows had undergone the label-free transform, the locked HC/MCI assessment
(21 HC and 20 MCI sample rows) gave AUC 0.576 (bootstrap 95% CI 0.387-0.776),
Brier 0.488 and ECE 0.488; target calibration was not fit. This weak and
uncertain ranking result with severe probability miscalibration was not treated
as positive replication or clinical validation. A separately locked,
development-only relative-expression top-scoring-pairs (TSP) transport audit
gave an internal AUC of 0.767 (95% CI 0.694-0.836), but AUCs of 0.548 (95% CI
0.486-0.612) in GSE140829 and 0.465 (95% CI 0.300-0.629) in GSE249477. Both
target domains were scored before endpoint assessment, including all 62
GSE249477 rows before diagnosis labels were parsed. Thus, relative-expression
ordering did not stably restore cross-platform transport.
Development-calibrated conformal sets achieved 0.669 mean external coverage at
nominal coverage 0.90 across 20 repeats. A phenotype-metadata-independent
domain audit, in which disease labels were excluded from preprocessing,
feature selection, thresholding and model fitting, separated development from
external expression domains with AUC 1.000 and permutation P = 0.00020.

These results define a practical reliability boundary for public blood-
transcriptome MCI claims. Same-family transfer and internal discrimination did
not imply independent transport, and computational uncertainty procedures did
not restore target-domain validity. The study therefore supports a leakage-aware
evaluation framework and a failure diagnosis, rather than a clinically deployable
MCI signature.

## Introduction

MCI occupies a clinically important but biologically heterogeneous position
between healthy ageing and dementia. Its outcomes vary from stable cognition to
progression toward Alzheimer disease (AD), creating a need for measurements that
can support risk stratification while respecting uncertainty in diagnosis and
trajectory. Blood-based molecular assays are attractive because they are more
accessible than brain tissue and may capture systemic or peripheral correlates of
neurodegeneration. However, the value of a blood molecular model depends not only
on its apparent discrimination, but also on whether the learned signal survives
new samples, new studies and new assay environments. The clinical and biomarker
background source remains pending final citation matching and is kept separate
from the methodological claims evaluated here.

Public transcriptome cohorts are useful for method development but create a
particularly difficult testing environment. Sample numbers are often modest,
the same study family may be represented by several accessions, expression
platforms differ in probe and gene coverage, and subject identifiers may be
incomplete or ambiguous. If feature selection, normalization, thresholding or
calibration uses information outside the training fold, an apparently high
internal estimate can include label information that would not be available at
deployment. Conversely, a rigorously fitted model can fail in an external cohort
because the expression domain changes. These are distinct problems and should be
measured separately, consistent with the distinction between internal,
internal-external and external validation and with standard definitions of
dataset shift (Steyerberg and Harrell, 2016; Moreno-Torres et al., 2012).

Clinical prediction reporting standards emphasize that model development,
validation, predictors, outcomes, performance and uncertainty should be
reported as separate elements. The same principle is relevant to molecular
biomarker studies, where a gene panel, an assay mapping and a disease endpoint
may each define a different estimand. TRIPOD, TRIPOD+AI, PROBAST and STROBE are
used here as reporting and bias-audit references, not as evidence that the
present model has passed a clinical validation threshold (TRIPOD, PMID
25560714; TRIPOD+AI, PMID 38626948; PROBAST, PMID 30596875; STROBE, PMID
18064739; all verified online on 2026-08-13).

Here we establish a leakage-aware public-data audit for blood-transcriptome MCI
models. We apply it to two AddNeuroMed development series, an independent GEO
series reserved for transport assessment, an RNA-seq stress-test cohort and a
separate progression branch. The primary question is not whether one model can
produce a high internal AUC. It is how much of the apparent signal remains after
fold-contained fitting, study-held-out transfer, exact-symbol mapping, clinical
baseline comparison, calibration, OOD assessment and phenotype-metadata-
independent domain diagnosis. This framing turns a negative external result
into an explicit reliability boundary while preserving historical computational
analyses as hypothesis-generating supplementary evidence.

## Results

### Cohort hierarchy and analysis units

The analysis was frozen before result generation and separated four evidence
roles. GSE63060 and GSE63061 formed the AddNeuroMed development family, with 155
and 142 blood bulk-microarray sample rows, respectively. Across the two series,
297 expression rows and 13,382 common genes were available for the domain audit;
the primary MCI-versus-HC endpoint contained 208 rows, including 69 MCI rows and
139 HC rows. The two accessions were treated as related study-family domains,
not as independent external populations.

GSE140829 contained 587 blood microarray rows, with explicit metadata labels of
249 HC, 134 MCI and 204 AD. The primary MCI-versus-HC analysis therefore used 383
rows. A validated sample-to-subject map was not established: candidate subject
prefixes repeated across samples and crossed diagnostic labels. They are
descriptive title fields, not subject identifiers, and must not be used for
cluster bootstrap, deduplication or subject-level confidence intervals. All
primary GSE140829 performance results are consequently reported as sample-row
transport audits, with subject independence unresolved.

GSE249477 was registered as a secondary independent cross-platform transport
audit, not as positive primary validation or replication. Its official series
labels comprised 21 HC, 20 MCI and 21 AD rows. GSE282742 remained a separate progression branch, comprising
116 sample rows, 86 unique subjects, 52 metadata-order representative P-MCI/S-MCI
records and 11 representative P-MCI-to-AD transitions. Public metadata did not
provide verified visit dates; the age/sex-adjusted representative-record analysis
therefore included 51 subjects (22 P-MCI and 29 S-MCI) and was not called a
baseline comparison. Because the progression estimand and primary
stage estimand differ, progression counts were not pooled with the primary
MCI-versus-HC analysis. GSE270454 was retained as an independent frozen-score
negative-boundary audit with 45 rows, including 25 HC, 10 MCI and 10 AD rows.
Before scoring, the count-column to GEO metadata join was verified as unique
for all 45 whole-blood RNA samples using descriptive cohort labels and terminal
sample ordinals.

The candidate compatibility audit also identified E-MTAB-6094 as unsuitable for
primary validation: its official BioStudies factor table contains 13 normal,
15 AD CDR1 and 7 AD CDR3 PBMC/mononuclear-cell samples, with no MCI group,
whereas the study description states 25 AD and 15 controls. It was excluded from
the primary endpoint and retained only as a compatibility-gate example.

An accession-identity audit prevented a further source of false replication.
The official GEO metadata cache contained 717 GSE63063 GSM records, exactly the
union of all 329 GSE63060 and all 388 GSE63061 GSM records, with no additional
GSE63063 record. GSE63063 was therefore retained only as a non-independent
aggregate audit and was excluded from model fitting, held-out testing and
independent-validation counts. It had no GSM overlap with GSE140829. We also
registered, but did not pool, three purified-neutrophil SubSeries
(GSE243743/GSE243744/GSE243745), the 60-row GSE208623 blood-leukocyte EPIC
methylation cohort, and two BCG-linked CSF/PBMC single-cell records
(GSE280720/GSE281606). GSE208623 supported a separate source-locked,
sex-adjusted MCI-versus-HC methylation association and descriptive gene-level
concordance with the frozen AddNeuroMed transcriptome reference. However, it
is a methylation assay rather than mRNA, age is unavailable and cell
composition is unadjusted; it therefore remains orthogonal contextual evidence
rather than validation of an unsorted whole-blood MCI-versus-HC model. The
other records' cell selection, longitudinal structure, exposure or stimulation,
and unresolved participant maps likewise preclude their use as validation.

### Leakage control separates modest internal optimism from transportability

The primary transcriptome comparator used exact-symbol expression intersections,
fold-contained feature screening and ridge logistic regression. In the primary
MCI-versus-HC development analysis, the nested internal AUC was 0.795 (95% CI
0.722-0.860). To quantify one plausible leakage route, we compared the clean
nested procedure with a repeated full-scope selection attack in which label-aware
feature selection was performed before diagnostic cross-validation. The attack
achieved mean AUC 0.822, whereas the clean nested diagnostic achieved 0.802,
yielding mean optimism of 0.0198 across repeated runs.

The optimism gap was measurable but modest relative to the subsequent external
failure. This distinction is important: fold-contained feature selection reduced
a known label-aware selection route, but it did not make the learned signal
transportable. The result supports a reliability audit in which leakage control
and external transport are reported as separate questions.

### Same-family transfer exceeds independent external transport

The transcriptome model transferred between the two AddNeuroMed accessions with
AUC 0.770 when trained on GSE63060 and tested on GSE63061 (n = 102), and AUC
0.730 in the reverse direction when trained on GSE63061 and tested on GSE63060
(n = 106). These estimates represent same-family study-held-out transfer, not
population-wide generalization.

The independent GSE140829 sample-row audit produced a primary transcriptome AUC
of 0.483 (95% CI 0.420-0.542; n = 383), below the no-discrimination reference
and not materially better than prespecified clinical baselines. Age-only,
sex-only and age-plus-sex baselines achieved AUCs of 0.518, 0.452 and 0.498,
respectively. Thus, the transcriptome model did not show incremental external
discrimination beyond these simple baselines in this audit. This pattern is
consistent with study-family dependence or expression-domain dependence, but it
does not identify the biological or technical cause.

### External uncertainty remains poor under domain displacement

The external transcriptome scores were poorly calibrated in the primary
sample-row audit. The Brier score was 0.344 and the ten-bin expected calibration
error was 0.340, compared with much lower descriptive calibration error in the
development assessment. Using a distance threshold defined only from development
data, 85.9% of external primary rows were above the training q95 OOD threshold.
The q90-q99 sensitivity analysis ranged from 57.4% to 96.3% for the mean
absolute standardized distance; training-only Ledoit-Wolf and robust-PCA
distances marked all external rows as OOD at each tested threshold. These values
are diagnostic quantities, not validated probabilities of domain membership.

We next tested two label-free procedures as stress tests. A within-sample
rank-then-z-score transform raised the external AUC to 0.605 (95% CI
0.538-0.663), but the Brier score remained 0.255 and the expected calibration
error remained 0.184. Development-only split-conformal calibration was also
evaluated at nominal coverage 0.90. Across 20 repeats, mean external coverage
was 0.669, with a repeated-run 2.5th-97.5th percentile interval of 0.650-0.779.
The nominal development guarantee therefore did not transfer to the external
sample-row audit under the observed shift.

The same separately locked transform was next assessed in GSE249477 as a secondary
independent cross-platform transport audit. We confirmed a one-to-one mapping
between all 62 series-title identifiers and TPM columns, then rank/z transformed
all 62 target rows before parsing target diagnosis labels. The final development
pipeline was fit only in AddNeuroMed; GSE249477 labels were used only once to
quantify the locked HC/MCI assessment after excluding 21 AD rows. In 41
assessment rows (21 HC, 20 MCI), AUC was 0.576 (bootstrap 95% CI 0.387-0.776),
Brier was 0.488 and ECE was 0.488. No target calibration, thresholding,
feature selection, fitting or transformation used target labels. This result
does not provide positive replication, clinical validation or evidence that the
rank/z transform solves domain shift.

We next performed a separately locked development-only relative-expression TSP
transport audit. In five-fold AddNeuroMed out-of-fold assessment, the method
gave an AUC of 0.767 (bootstrap 95% CI 0.694-0.836). Ten non-overlapping,
MCI-oriented gene pairs were then selected from all 208 AddNeuroMed HC/MCI rows
and frozen for target scoring. All 587 GSE140829 expression rows were scored
before endpoint filtering; the final 383-row HC/MCI sample-row assessment gave
an AUC of 0.548 (95% CI 0.486-0.612). All 62 GSE249477 rows were likewise
scored before diagnosis labels were parsed; its locked 41-row HC/MCI assessment
gave an AUC of 0.465 (95% CI 0.300-0.629). GSE249477 target calibration was not
fit. The relative-expression audit therefore did not provide conclusive
transport in GSE140829 or reproduce discrimination in the secondary RNA-seq
assessment; it is not external validation or evidence of clinical deployment.

### Direction-concordant signals retain development-family consistency but not external transport

We next tested whether an explicitly cross-study stable feature construction
would change this boundary. Before target diagnosis metadata were accessed, we
restricted the candidate universe to 11,575 exact symbols available in both
AddNeuroMed studies and both target platforms. In synchronized five-fold outer
partitions, standardized MCI-minus-HC effects were estimated independently in
the GSE63060 and GSE63061 training partitions; only non-zero effects with the
same direction were eligible. The top 100 genes were ranked by their minimum
absolute study effect, with deterministic ties, and a fixed ridge model was fit
only to the combined paired training rows.

This procedure produced AUC 0.840 (95% CI 0.753-0.913) in the GSE63060 outer
folds and 0.779 (95% CI 0.677-0.870) in GSE63061. The pooled development OOF
AUC was 0.819; a 100-repeat label-permutation path that repeated feature
selection and outer-fold fitting gave P=0.0099. These values support
development-family directional consistency, but do not establish a new
independent validation domain.

After a single final selection on all AddNeuroMed MCI/HC rows, all 587
GSE140829 and all 62 GSE249477 expression rows were scored before their disease
labels were used for assessment. The locked MCI-HC AUC was 0.541 (95% CI
0.479-0.598; n=383) in GSE140829 and 0.524 (95% CI 0.500-0.575; n=41) in
GSE249477. External labels were not used for mapping, feature selection,
preprocessing, fitting, prediction calibration, thresholding or adaptation; as
with the other audits, they were used after prediction lock only for descriptive
performance and calibration assessment. Direction-concordant selection therefore
did not restore external transport or establish a clinically valid biomarker.

### Gene-level effects are related-study concordant but cross-platform discordant

To distinguish model failure from instability in the underlying gene-level
contrast, we performed a non-predictive effect-concordance audit in the same
11,575-symbol universe. The universe was fixed from expression mapping before
external diagnosis labels were accessed. Independently estimated MCI-minus-HC
Hedges' g values were substantially concordant between GSE63060 and GSE63061
(Spearman correlation 0.673; 71.8% nonzero sign concordance). The unweighted
AddNeuroMed descriptive reference remained positively correlated with GSE140829
(Spearman 0.486; 66.2% sign concordance), but GSE140829 retained its sample-row
unit and unresolved subject independence.

In contrast, the GSE249477 comparisons were directionally discordant: Spearman
correlations with GSE63060, GSE63061 and GSE140829 were -0.182, -0.174 and
-0.191, respectively, with sign concordance of 43.8%, 44.4% and 42.6%. These
comparisons used 10,923 estimable genes because 652 shared genes lacked the
within-group information needed for a valid standardized effect in GSE249477.
We reported Cochran Q/I-squared only for the two AddNeuroMed series as a
descriptive within-family statistic; no four-cohort pooled effect or
heterogeneity P value was calculated. Thus, the audit quantifies a marked
cross-platform reliability boundary without attributing it to a biological,
cellular or technical cause (Fig. 7).

### The frozen 12-gene comparator is bounded by exact mapping

The historical frozen comparator contained GPI, PHF15, VAT1, EIF3B, UCP2,
LAT1-3TM, C22ORF39, USP38, ZSWIM6, PRDM10, RP2 and GCA. The exact symbols were
available in the development data, but only 10 of 12 mapped exactly in GSE140829
and GSE249477. The absent symbols differed by cohort, and no alias substitution
was permitted. The full score was therefore blocked externally; partial scores
were retained only as sensitivity analyses.

The partial frozen score gave AUC 0.580 for MCI versus HC in GSE140829 and AUC
0.410 in the GSE249477 RNA-seq stress test. In the independent GSE270454
boundary audit, only 10 of 12 symbols mapped, the verified row-level
MCI-versus-HC AUC was 0.404 (bootstrap 95% CI 0.204-0.612), and the
descriptive age- and sex-adjusted MCI-versus-HC estimate was -0.147
(P = 0.373; n = 45). These results do not constitute complete-panel replication.
They establish that exact-symbol coverage and score provenance are part of the
validity claim, rather than implementation details that can be silently repaired
after seeing external results.

### Domain separation is detectable without disease labels in model fitting

To test whether the external failure was accompanied by a measurable expression
domain boundary, we trained a domain-only classifier to distinguish the 297
AddNeuroMed development rows from the 587 GSE140829 external rows. Disease
metadata were not used for preprocessing, feature selection, threshold choice or
domain-model fitting. They were used only for row alignment and for post hoc
phenotype-stratified sensitivities. In a repeated 10-by-5-fold audit over 13,382
common genes, domain AUC and balanced accuracy were both 1.000, with permutation
P = 0.00020.

Post hoc HC, MCI and AD strata each retained domain AUC 1.000. This strengthens
the interpretation that the two expression domains are strongly separable in
this dataset, but it does not show that the domain classifier identified a
disease mechanism or prove why the MCI model failed. The finding is best treated
as evidence of severe domain displacement that must be addressed before clinical
transportability can be assessed.

## Discussion

The central result is a separation between apparent discovery performance and
transportability. A fold-contained transcriptome model discriminated MCI from HC
within the AddNeuroMed development family, and the two related accessions
supported study-held-out AUCs between 0.730 and 0.770. Yet the same locked
external assessment failed to discriminate MCI from HC in GSE140829, with poor
calibration and a high training-derived OOD rate. This is not a contradiction.
It is the expected pattern when a model captures study-family expression
structure, disease-correlated technical structure, or a mixture of both.

The leakage attack provides a concrete warning but not the whole explanation.
The average optimism gap of 0.0198 shows that one label-aware selection route can
inflate internal estimates. The much larger difference between same-family
transfer and independent transport indicates that leakage control alone is
insufficient. A model can be statistically clean within a development protocol
and still fail when the expression domain changes.

The domain-only result provides a measurable diagnosis of that boundary. Perfect
domain separation persisted after disease labels were removed from model fitting,
and it was observed within post hoc disease strata. This argues against an
interpretation based only on gross case-control imbalance, but several rival
explanations remain possible, including platform processing, cohort recruitment,
batch structure, medication, demographic composition, RNA quality and
unresolved repeated sampling. The current evidence cannot distinguish these
causes. The correct claim is therefore domain displacement, not a biological
explanation of model failure.

The sensitivity analyses illustrate why domain adaptation should not be judged by
one discrimination number. A sample-wise rank/z-score transformation partially
recovered ranking, but calibration remained limited. A development-only
relative-expression TSP audit also retained internal discrimination but did not
stably transport across the microarray and RNA-seq assessments. Conformal sets
calibrated in the development domain did not achieve nominal external coverage.
These analyses are useful failure diagnostics under shift; they do not create a
clinical guarantee. The method sources are recorded in the verified register:
pairwise mRNA comparison and relative-expression classifier methods (PMID
16646797; PMID 23496976), dataset shift (Moreno-Torres et al., 2012) and
distribution-free predictive inference (Lei et al., 2018).

The frozen 12-gene panel also demonstrates a provenance boundary. Exact mapping
was complete in development but incomplete in each external assay, and the
historical source record does not yet provide a complete publication-level
discovery provenance record. The panel can therefore be retained as a fixed
audit comparator, but it cannot be described as a validated independent
signature. Likewise, MR, virtual knockout, LINCS perturbation, network
pharmacology, docking, molecular dynamics, ADME and cell-level simulations remain
computational hypothesis generation. They may prioritize future experiments, but
they cannot rescue an expression-domain failure or establish causality, binding,
efficacy or clinical utility.

### Limitations

First, the GSE140829 subject-to-sample map remains unresolved. The primary
external confidence intervals and error rates are therefore sample-row audits;
they may not represent independent-subject uncertainty if repeated rows are
present. Second, the independent cohorts do not provide a complete exact-symbol
mapping for the frozen panel, and the RNA-seq cohort was analysed only as a
cross-platform stress test. Third, the development sample size is modest and the
observed domain classifier may reflect multiple correlated acquisition or cohort
factors. Fourth, no prospective cohort, clinical utility analysis, decision-curve
assessment or treatment-impact study was performed. Fifth, the progression branch
has a distinct estimand, no verified visit dates and no FDR-positive gene-level
associations. Its supplementary dual-source-mapped Reactome analysis yielded
exploratory competitive pathway associations, but these do not identify a cell
composition change, establish a mechanism or provide independent predictive
validation of the primary MCI-versus-HC model. Finally, online verification of the background and
method references remains incomplete in the current environment and must be
resolved before submission.

### Implications and next study

The immediate implication is methodological. Public blood-transcriptome MCI
studies should report the unit of analysis, exact assay mapping, fold-contained
selection, study-family structure, clinical baselines, calibration, OOD
diagnostics, uncertainty assumptions and independent transport results in one
evidence chain. For this model family, a credible next step is not further
post-hoc tuning of the current public cohorts. It is a prespecified, multi-site
subject-level cohort with a locked assay mapping, explicit repeated-sample
handling, blinded external scoring, domain-aware calibration and clinical utility
evaluation. The present package provides the audit criteria and failure modes
that such a study should resolve.

## Conclusion

Strict leakage control improved the interpretability of public blood-transcriptome
MCI modelling, but it did not convert internal or same-family discrimination into
independent transportability. In the present audit, the external expression
domain was strongly separable from development data, external calibration and
conformal coverage were poor, and exact frozen-panel mapping was incomplete.
The defensible contribution is a reproducible reliability boundary for MCI
biomarker claims and a set of auditable failure diagnostics. The evidence does
not support a clinically validated diagnostic signature.

## Methods

### Study design and data sources

This retrospective public-data study used GEO expression matrices and official
sample metadata for GSE63060, GSE63061, GSE140829 and GSE249477. The GSE140829
GEO title exactly matched the registered bioRxiv preprint
`10.1101/2019.12.13.875112`; this provenance link does not resolve the
sample-to-subject map or change its sample-row transport-audit role. GSE282742 was
retained as a separate progression branch, and GSE270454 as a frozen-score
boundary audit. GSE270454 was verified against the official NCBI GEO Entrez
summary (GSE270454, GPL24676) and the linked publication ([PMID 40779089](https://pubmed.ncbi.nlm.nih.gov/40779089/));
all main analysis accessions have an accession-level source audit before
release. Official accession records and cohort-paper or preprint provenance
status are listed in the V5 source-verification register.

Candidate cohorts were screened by assay, specimen, label, mapping and analysis
unit before inclusion. Direct blood bulk expression, PBMC expression, induced
cell-response assays, non-coding RNA data and progression data were registered as
separate evidence domains. Cohorts were not pooled merely to increase sample
count. GSE18309 was excluded from primary validation because its compatible
cohort was too small. GSE312139 was retained as orthogonal induced-neuronal
response context rather than direct blood expression. E-MTAB-6094 was excluded
from the primary endpoint because it is PBMC/mononuclear-cell AD-versus-control
data with no MCI group and a factor-table versus description count conflict.
GSE63063 was accession-audited as a complete GSE63060/GSE63061 aggregate and
excluded from all independent-validation and model-input roles. Candidate
neutrophil SubSeries and BCG-linked CSF/PBMC single-cell records were kept only
as explicitly non-pooled supplementary contexts.

An explicit official GEO metadata query for `"mild cognitive impairment"` and
`blood` was saved on 2026-08-13. It returned 38 unique GSE records, including
18 that were not already registered. No newly unregistered record met all
prespecified requirements for a human, direct peripheral whole-blood bulk-mRNA,
independent MCI-versus-HC comparator. This is a saved query scope, not a claim
that all public archives have been exhausted; the exact artifacts and screening
record are retained in the candidate-discovery note.

### Endpoints and analysis units

The primary endpoint was MCI versus HC. MCI versus AD was a key secondary
endpoint. AD versus HC was retained as contextual analysis. Progression was
defined separately as the P-MCI versus S-MCI contrast in GSE282742 and was not
used to validate the primary endpoint. Sample rows were the analysis unit for
the development and expression-domain audits. GSE140829 results were explicitly
labelled sample-row audits because subject independence was unresolved. The
progression branch used one metadata-order representative record per unique
subject because verified visit dates were absent, and paired representative
subject transitions where applicable.

### Leakage-free transcriptome comparator

Expression features were restricted to exact-symbol intersections after
source-platform mapping. Probe aggregation and other unsupervised mapping steps
were performed without target labels. Within every development fold, feature
screening, imputation, scaling and model fitting were fitted using training rows
only. Feature screening used `SelectKBest(f_classif, k=min(200, n_features))`.
The classifier was ridge logistic regression with fixed `C = 0.2`. The final
external model was fitted once on all eligible development rows after internal
nested evaluation was completed. External labels were not used for fitting,
feature selection, preprocessing, calibration, threshold selection, OOD
threshold choice or adaptation.

Age-only, sex-only and age-plus-sex logistic baselines were prespecified and
reported alongside the transcriptome comparator. Calibration intercept and slope
were descriptive evaluation-set assessments fitted after predictions were locked;
they were not used to recalibrate external scores. AUC uncertainty was estimated
with bootstrap intervals over the stated analysis unit. ECE used fixed
probability bins, with the number of non-empty bins recorded in the result table.

### Leakage attack and study-held-out transfer

The leakage attack reproduced a plausible full-development, label-aware feature
selection route before repeated diagnostic cross-validation. It was compared
with the fold-contained nested procedure using the same primary endpoint and
model family. Optimism was defined as the attack AUC minus the clean nested AUC
within each repeat, then summarized across repeats. Study-held-out transfer was
performed bidirectionally between GSE63060 and GSE63061. Because both accessions
belong to AddNeuroMed, these results were interpreted as same-family transfer,
not independent population validation.

### External transport assessment

GSE140829 was held out from all model-fitting decisions. Its official series
matrix and GPL15988 platform mapping were used for expression alignment. The
primary endpoint was restricted to rows explicitly labelled MCI or HC. Exact
sample-to-subject independence was not assumed. External labels were used only
for post hoc performance, calibration, error and sensitivity assessment.

GSE249477 was processed as a secondary independent cross-platform transport
audit. It used blood bulk RNA-seq and was not combined with the primary
microarray estimand. After a strict 62/62 title-to-TPM identifier check, all 62
target rows underwent within-sample rank then z-score transformation before
target diagnosis labels were parsed. The AddNeuroMed-only pipeline retained
median imputation, fold-contained feature selection, scaling and ridge logistic
regression. GSE249477 labels were used only for the one locked 21-HC versus
20-MCI assessment after 21 AD rows were excluded. They were not used for target
transformation, model fitting, feature selection, thresholding or calibration;
no target calibration was fit.

### Development-only relative-expression TSP transport audit

As a separately locked methodological extension, we tested whether a
relative-expression classifier could reduce platform-scale dependence without
using target-domain fitting. Pairwise mRNA comparison and related
relative-expression classifiers motivate this method family (PMID 16646797;
PMID 23496976), but these records do not provide evidence for MCI clinical
utility. Within each of five AddNeuroMed training folds, `f_classif` ranked
development-only genes, the top 100 finite-score genes were retained with
deterministic tie handling, and the ten largest MCI-versus-HC pair-ordering
gaps were greedily selected without reusing a gene. Each pair was oriented so
that its MCI-enriched ordering contributed one vote; a sample score was the
mean oriented vote across the frozen pairs. Final pairs were selected only in
the complete 208-row AddNeuroMed HC/MCI development set after out-of-fold
assessment.

No target-domain normalization, pair selection, classifier fitting, threshold,
calibration or adaptation was fitted. GSE140829 was scored across all 587
retained rows before the HC/MCI endpoint was used for the final 383-row
assessment. GSE249477 was scored across all 62 expression rows before diagnosis
labels were parsed, then the frozen 41-row HC/MCI assessment was calculated.
Target labels were used only after predictions were fixed for aggregate endpoint
metrics. Calibration summaries for GSE140829 were descriptive assessments after
prediction locking; calibration intercept and slope were intentionally not fit
for GSE249477 target scores.

### OOD, selective prediction and conformal stress tests

The primary OOD distance was derived from training-only standardized expression
features. Thresholds were fixed at development q90, q95 and q99 values. A
training-only Ledoit-Wolf Mahalanobis distance and a robust-PCA distance were
added as sensitivity definitions. External labels were used only to assess
errors and discrimination after the distances and thresholds were locked.

Selective-risk curves used confidence thresholds indexed only by development
predictions. Split-conformal prediction used a development-only calibration
split with nominal coverage 0.90 and 20 repeated fits. External labels were used
only to evaluate coverage, set size, singleton error and AUC. No claim of
distribution-free target-domain coverage was made under the observed shift.

### Domain-separability audit

The domain classifier used a binary domain label, with AddNeuroMed development
rows as one class and GSE140829 external rows as the other. All expression rows
were retained for the primary domain audit. Phenotype metadata were not used for
preprocessing, feature selection, thresholding or model fitting. They were used
for row alignment and post hoc HC, MCI and AD stratum definitions. Domain
classification used repeated 10-by-5-fold cross-validation, 200 training-fold
selected features and a permutation test with 5,000 permutations. The audit was
interpreted as evidence of expression-domain separation, not as a disease model
or biological causal analysis.

### Frozen-panel provenance and mapping

The historical 12-gene comparator was frozen before V5 evaluation. Its exact
symbols and fixed direction-based audit weights were preserved without
reselection, coefficient refitting, alias substitution or replacement of
`LAT1-3TM` by `SLC7A5`. If an assay lacked an exact symbol, the complete score was
blocked and a partial score was reported only as a sensitivity audit. The
provenance ledger records the historical source table and fixed weights while
marking publication-level discovery provenance as incomplete.

### Reproducibility and statistical reporting

All result tables, figure source data, scripts and aggregate hashes are stored
in the isolated V5 package. The package does not copy raw matrices or sample
identifiers into the release directory. Six main figures are exported as PDF,
SVG and 600-dpi PNG variants, with a figure manifest and source-data directory.
The internal QA contract checks script compilation, result dimensions, endpoint
boundaries, mapping rules, analysis-unit labels, figure resolution, blank-image
checks, source-data presence and manifest hashes.

## Data availability

The primary and supplementary analyses use publicly available GEO accessions
[GSE63060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060),
[GSE63061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061),
[GSE140829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829),
[GSE249477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249477),
[GSE282742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282742) and
[GSE270454](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270454).
Official source URLs, platform records, original-study links and verification
dates are recorded in `source_verification_register.md`. The V5 package
contains aggregate ledgers, anonymized prediction rows, source-data tables,
scripts and file hashes; it does not redistribute raw expression matrices or
sample identifiers. `[AUTHOR INPUT NEEDED: final public repository DOI or
archive URL and release license statement]`

## Code availability

The reproducible V5 analysis package contains the protocol, analysis scripts,
result tables, figure-generation script, figure source data, manifest hashes and
QA tests. The direct Python dependency pins and verified Python 3.10 execution
record are provided in `requirements-v5.txt` and
`docs/reproducibility_environment_v1.md`; `99_QA/check_runtime_environment.py`
checks those direct versions before the full QA. `[AUTHOR INPUT NEEDED: public
repository URL, immutable commit/tag, resolved target-platform environment
export and archival DOI]`

## Ethics statement

This study used previously released, de-identified public expression matrices
and metadata. No new participants were recruited and no new identifiable data
were accessed. The ethics and consent statements of each source study should be
checked against the official GEO record and original publication before release.

## Figure legends

### Fig. 1 | Design and cohort hierarchy

The frozen endpoint and model definitions are applied to the AddNeuroMed
development family, followed by the GSE140829 independent external transport
assessment and the GSE249477 secondary independent cross-platform transport
audit. GSE282742 is a separate
progression branch and GSE270454 is a frozen-score boundary audit. Counts are
sample rows unless a unique-subject count is explicitly shown. The primary
endpoint is MCI versus HC. This workflow panel provides study-design context and
does not substitute for the result panels.

### Fig. 2 | Leakage and transportability

(a) Repeated five-fold diagnostic comparison of fold-contained selection with a
full-development label-aware selection attack. (b) Bidirectional transfer
between the two AddNeuroMed series and transport assessment in GSE140829. (c)
Primary MCI-versus-HC transcriptome AUC in development and GSE140829. Points show
AUC and bars show bootstrap 95% intervals. The external result is a sample-row
transport audit, not positive clinical or subject-level validation.

### Fig. 3 | Primary-endpoint external uncertainty

All panels are restricted to GSE140829 MCI-versus-HC rows (n = 383; HC = 249,
MCI = 134). (a) Five-bin descriptive calibration curves for the age-plus-sex and
transcriptome models. (b) Training-derived standardized distance by label with
the development q95 threshold. This distance is diagnostic, not a validated
domain-membership probability. (c) Selective error against external coverage
using thresholds indexed only by development confidence quantiles. External
labels are assessment-only.

### Fig. 4 | Exact-symbol mapping boundary

The frozen signed 12-gene score maps exactly to 12/12 development symbols but
only 10/12 symbols in GSE140829, GSE249477 and GSE270454. Missing symbols are
not substituted. The GSE270454 count-to-GEO join was verified for 45/45
whole-blood RNA samples before scoring. Partial external scores are sensitivity
or boundary analyses. AUC intervals are bootstrap 95% intervals over sample
rows.

### Fig. 5 | Development-family signal exceeds a label null but does not restore external transport

(a) The label-free within-sample rank/z-score transform gives AUC 0.605 in the
GSE140829 sample-row audit and AUC 0.576 (bootstrap 95% CI 0.387-0.776; n=41)
in the secondary independent cross-platform GSE249477 transport audit. The
development-only relative-expression TSP transport audit gives AUC 0.548
(95% CI 0.486-0.612) in GSE140829 and 0.465 (95% CI 0.300-0.629) in GSE249477.
TSP pairs were selected solely in AddNeuroMed; all 587 GSE140829 rows and all
62 GSE249477 rows were scored before endpoint assessment, and GSE249477 labels
were parsed only afterward. (b) GSE140829 AUC, ECE and Brier are shown on
metric-specific axes and are not a common clinical utility scale. (c)
Development-calibrated split-conformal prediction at nominal coverage 0.90 is
evaluated in the GSE140829 external sample-row audit; mean coverage is 0.669
across 20 repeats. These panels are stress tests, not external validation,
clinical deployment evidence, or evidence that rank/z or TSP solves domain
shift.

### Fig. 6 | Phenotype-metadata-independent domain separability

The domain classifier distinguishes AddNeuroMed development rows from GSE140829
external rows using expression values and a binary domain label. Phenotype
metadata were discarded before domain preprocessing, feature selection and model
fitting, but were used for row alignment and post hoc stratum definition. (a)
Receiver-operating-characteristic curve. (b) Predicted probability distributions
by domain. (c) AUC across ten repeated five-fold cross-validation runs. Aggregate
domain AUC is 1.000, with permutation P = 0.00020. This is evidence of domain
separation, not a disease classifier, and GSE140829 remains a sample-row audit.

### Fig. 7 | Cross-cohort MCI-HC effect concordance and descriptive heterogeneity

(a) Pairwise Spearman correlations of independently estimated cohort-wise
MCI-minus-HC Hedges' g in a 11,575-symbol universe locked before external
diagnosis-label access. (b) Pairwise nonzero sign concordance; GSE249477 pairs
use 10,923 estimable genes. (c) Descriptive Cochran Q/I-squared distribution
for GSE63060/GSE63061 only. GSE140829 is excluded from pooled inference because
it is a sample-row audit with unresolved subject independence. No pooled
four-cohort effect or heterogeneity P value is reported.

## Supplementary analysis index

The historical cis-eQTL-to-AD MR and colocalization output is retained only as
an exclusion-audit case: `audit_legacy_mr_method_boundary.py` records failures
of the V5 effect-scale, LD-clumping, gene-locus, phenotype and rerun-provenance
gates. It is not causal MCI evidence. Other historical virtual knockout, LINCS,
network pharmacology, docking, molecular dynamics, ADME, scGen/Augur and
cell-level simulation analyses remain in the source project and prior
integration package. They should be reported only as `computational hypothesis
generation` or `orthogonal context`. They cannot be used to claim causal MCI
effects, experimental mechanism, therapeutic efficacy, direct blood-mRNA
validation or clinical utility. The supplementary package should link every
panel to its own source table, script, software version, random seed and
evidence boundary.

## Reference and verification register

The following source-verification register is current as of 2026-08-14. It
distinguishes verified method/reporting records, cohort-paper records and the
title-matched GSE140829 preprint from remaining frozen-panel provenance gaps.

1. Frozen 12-gene panel discovery provenance: `PENDING_FINAL_MATCH`.
   Main analysis accessions have accession-level source audits; GSE140829 is
   linked to a title-matched preprint, not represented as peer-reviewed or
   subject-level validation evidence.
2. Collins et al., TRIPOD, PMID 25560714, DOI `10.7326/M14-0697`:
   `VERIFIED_ONLINE`.
3. Collins et al., TRIPOD+AI, PMID 38626948, DOI
   `10.1136/bmj-2023-078378`: `VERIFIED_ONLINE`.
4. Wolff et al., PROBAST, PMID 30596875, DOI `10.7326/M18-1376`:
   `VERIFIED_ONLINE`.
5. von Elm et al., STROBE, PMID 18064739, DOI
   `10.1016/S0140-6736(07)61602-X`: `VERIFIED_ONLINE`.
6. Steyerberg and Harrell, internal, internal-external and external validation,
   DOI `10.1016/j.jclinepi.2015.04.005`: `VERIFIED_ONLINE`.
7. Moreno-Torres et al., dataset shift, DOI `10.1016/j.patcog.2011.06.019`:
   `VERIFIED_ONLINE`.
8. Lei et al., distribution-free predictive inference, DOI
   `10.1080/01621459.2017.1307116`: `VERIFIED_ONLINE`.
9. Pairwise mRNA comparison classifier, PMID 16646797:
   `VERIFIED_ONLINE`; method-family rationale only.
10. AUREA relative-expression classifier methods, PMID 23496976:
    `VERIFIED_ONLINE`; method-family rationale only.
11. Bootstrap/ECE implementation details, frozen-panel discovery provenance and
   retained supplementary methods: `PENDING_METHOD_AND_PANEL_SOURCE_MATCH`.

## Internal evidence binding

The current numerical claims are traceable as follows:

| Claim family | V5 source | Figure / section |
|---|---|---|
| Development, external and analysis-unit counts | `01_DATA_AUDIT/data_ledger.csv`; `00_PROTOCOL/domain_registry.csv` | Fig. 1; Results cohort section |
| Nested AUC, baselines, external calibration and OOD | `03_RESULTS/metrics.csv` | Fig. 2-3; Results transport sections |
| Leakage optimism | `03_RESULTS/leakage_attack_summary.csv` | Fig. 2; Results leakage section |
| Same-family transfer | `03_RESULTS/study_heldout_metrics.csv` | Fig. 2; Results transfer section |
| Rank transport, relative-expression TSP transport, and conformal stress tests | `03_RESULTS/domain_invariant_sensitivity.csv`; `gse249477_rank_transport_metrics.csv`; `gse249477_rank_transport_audit.json`; `tsp_transport_metrics.csv`; `tsp_transport_selected_pairs.csv`; `tsp_transport_audit.json`; `split_conformal_summary.csv` | Fig. 5; Table S13 |
| OOD distance sensitivity | `03_RESULTS/ood_distance_sensitivity.csv`; `ood_distance_audit.json` | Fig. 3; Methods OOD section |
| Frozen panel and mapping boundary | `03_RESULTS/panel_provenance_ledger.csv`; `metrics.csv`; `gse249477_stress_test.json`; `gse270454_join_audit.json`; `gse270454_boundary_audit.json` | Fig. 4 |
| Domain separability | `03_RESULTS/domain_separability_summary.csv`; `domain_separability_strata.csv`; `domain_separability_audit.json` | Fig. 6 |
| Cross-cohort effect concordance | `03_RESULTS/cross_cohort_effect_pairwise_summary.csv`; `cross_cohort_effect_heterogeneity.csv`; `cross_cohort_effect_audit.json` | Fig. 7; Table S15 |
| Claim restrictions | `05_MANUSCRIPT/claim_evidence_ledger.md` | Entire manuscript |
