# MCI Public Reliability Analysis V5

This is an isolated, public-data-only analysis package created on 2026-08-12.
It does not overwrite V4, the frozen `MCI_DOMAIN_STABLE_20260811` audit package,
or any canonical MCI input root. V4 remains frozen as the preceding release.

The package has two deliberately separate layers:

1. An exact-symbol frozen 12-gene comparator. It is blocked for the complete
   GSE140829 external cohort because two exact symbols are absent from the
   official GPL15988 mapping.
2. A leakage-free transcriptome comparator. Its features are selected only
   inside development folds, then the final model is evaluated once on the
   GSE140829 external transport-assessment cohort.

The package also contains bidirectional study-held-out transport between the
two development series, a repeated diagnostic leakage attack, a corrected
public-data ledger separating sample rows from unique subjects, seven main
figures plus one supplementary figure with source data and hashes, and a
claim-evidence ledger that blocks
unsupported clinical and causal wording.

V5 adds a phenotype-blind domain-separability audit and a per-gene frozen-panel
provenance ledger. The domain classifier uses a binary domain label
(AddNeuroMed development versus GSE140829 external); phenotype metadata are
used only for public-expression row alignment and are discarded before domain
preprocessing, feature selection and model fitting. The provenance ledger records the historical source table and
fixed direction-based audit weights without presenting them as a complete
independent discovery record.

This version adds a candidate compatibility registry for newly discovered GEO
series. Direct blood transcriptome data, PBMC data, cell-response assays,
non-coding RNA data, and progression data are registered separately and are
not pooled across incompatible estimands. GSE18309 is excluded from primary
validation because it contains only three MCI and three control samples; the
new GSE312139 iCAP series is retained as orthogonal context because its readout
is an induced neuronal response to plasma, not direct blood expression.
E-MTAB-6094 is excluded from the primary endpoint because it is PBMC/mononuclear-
cell AD-versus-control data with no MCI group, and its official factor table
conflicts with the study-description sample totals. The field-level audit is
retained as a compatibility-gate example.

The official GEO metadata cache also shows that GSE63063 is a complete GSM
aggregate of the AddNeuroMed development records: all 329 GSE63060 and all 388
GSE63061 GSM accessions occur in GSE63063, with no extra GSE63063 GSM records.
It is excluded from independent validation and model inputs. The audit also
records zero GSM overlap between GSE63063 and GSE140829. Candidate neutrophil
SubSeries GSE243743, GSE243744 and GSE243745 are separately registered as a
single-study-family, cell-type-specific context only; the BCG-linked
GSE280720/GSE281606 records are longitudinal PBMC/CSF single-cell intervention
context only. None is pooled into the primary whole-blood bulk MCI-HC endpoint.

GSE208623 is separately audited as an independent-source Chinese peripheral-
blood-leukocyte EPIC DNA-methylation cohort (20 HC, 20 MCI and 20 AD). It is
excluded from all mRNA model fitting, tuning, calibration, validation and
pooling. The source-locked processed beta matrix and EPIC annotation support a
separate sex-adjusted MCI-HC CpG association in 20 MCI and 20 HC rows, with AD
excluded before modeling. Its gene-level signed-Z aggregates have descriptive
concordance with the frozen AddNeuroMed transcriptome-effect reference; this is
orthogonal association context, not external mRNA validation, replication,
methylation-expression mediation, cell-composition proof or causal evidence.
Age and blood-cell composition are not adjusted.

It also includes label-free sample-wise rank/z-score normalization and
development-only split-conformal uncertainty sensitivities. The rank/z
sensitivity is additionally assessed in GSE249477 as a secondary
independent cross-platform transport audit: 41 HC/MCI assessment rows yielded
AUC 0.576 (bootstrap 95% CI 0.387-0.776), Brier 0.488 and ECE 0.488, with no
target calibration. All 62 target rows were transformed before diagnosis labels
were parsed; labels were used only for the locked final assessment. These are
stress tests, not target-domain adaptation or positive replication. A separate OOD distance
audit compares training-only standardized, Ledoit-Wolf Mahalanobis, and robust
PCA distances at development q90/q95/q99 thresholds; external labels are used
only for post hoc error assessment.

This version also adds a separately locked, development-only relative-expression
top-scoring-pairs (TSP) transport audit. Ten disjoint, MCI-oriented gene pairs
are selected only in AddNeuroMed, then scored without target fitting,
calibration, threshold selection or adaptation. The out-of-fold development AUC
is 0.767 (95% CI 0.694-0.836); target AUCs are 0.548 (95% CI 0.486-0.612) in
GSE140829 and 0.465 (95% CI 0.300-0.629) in GSE249477. All 587 GSE140829 rows
and all 62 GSE249477 rows are scored before endpoint assessment. TSP is a
method stress test, not external validation, clinical deployment, or a solution
to domain shift.

V5 also includes a cross-study stable disease-signal audit. A mapping-only
candidate universe of 11,575 exact symbols common to both development studies
and the two target platforms was declared before target diagnosis access. In
synchronized GSE63060/GSE63061 outer folds, only genes with direction-concordant
MCI-minus-HC effects in both training partitions were eligible for a fixed
100-gene ridge model. The development-family stable model scored AUC 0.840 in
GSE63060 and 0.779 in GSE63061; a complete 100-repeat label-permutation path
gave p=0.0099 for pooled development OOF AUC 0.819. However, locked transport
was AUC 0.541 (95% CI 0.479-0.598) in GSE140829 and 0.524 (95% CI 0.500-0.575)
in GSE249477. This documents development-family consistency, not independent
validation or recovery of external transport.

The primary endpoint is MCI versus healthy control. MCI versus AD is a key
secondary endpoint. Progression remains a separate endpoint.

Run the analysis with the bundled workspace Python runtime:

`E:\AI_Research_Assistant\scripts\utils\python.cmd 02_ANALYSIS\run_pairwise_external.py`

Before reproducing on another machine, use Python 3.10 and install the direct
dependency pins in `requirements-v5.txt`. The observed V5 runtime and the
remaining archive-release requirement are documented in
`docs/reproducibility_environment_v1.md`.

Run `99_QA\check_runtime_environment.py` before the full QA to verify that the
direct package versions match the V5 execution record.

Then run the V5 additions:

```text
02_ANALYSIS\build_panel_provenance.py
02_ANALYSIS\audit_gse270454_join.py
02_ANALYSIS\run_gse270454_locked12_joined.py
02_ANALYSIS\run_domain_separability_audit.py
02_ANALYSIS\run_tsp_relative_expression_transport.py
02_ANALYSIS\run_cross_study_stable_signal.py
02_ANALYSIS\audit_emtab6094_compatibility.py
02_ANALYSIS\make_figure_pack.py
99_QA\run_full_qa.py
```

The script writes aggregate results and hashes only under this package. It
does not copy raw data or sample identifiers into the package.

## Current evidence boundary

The primary MCI-versus-HC transcriptome model has nested internal AUC about
0.80 in the AddNeuroMed development family, but AUC about 0.48 in the
GSE140829 external transport assessment with a high OOD rate. The GSE140829
result is a sample-row audit because a validated sample-to-subject map is not
available. GSE249477 is a small RNA-seq secondary cross-platform transport
audit with incomplete exact-symbol mapping. Its rank/z result is weak and
uncertain and has severe probability miscalibration, and the development-only
relative-expression TSP audit has AUC 0.465 (95% CI 0.300-0.629). Neither is
external validation. These results
support a reliability/transportability manuscript, not a validated clinical
biomarker claim.

The progression branch uses GSE282742 as a separate exploratory association:
116 sample rows, 86 unique subjects, 52 metadata-order representative P-MCI/S-MCI
records (22 P-MCI, 30 S-MCI), and 51 age/sex-complete records (22 P-MCI, 29
S-MCI). Public metadata does not provide verified visit dates, so these are not
called baseline samples. Eleven representative P-MCI-to-AD pairs provide only a
paired-transition sensitivity; neither analysis has FDR-positive genes. The
the pathway layer uses a version-pinned Reactome gene-set registry and remains
an exploratory association summary. This is not primary-stage validation or
prediction evidence.

GSE270454 is an independent small whole-blood RNA-seq boundary audit. Its
count-column to GEO metadata join was verified as unique for 45/45 RNA samples
using descriptive cohort labels and terminal sample ordinals. The frozen panel
mapped 10/12 exact symbols, so it remains a partial frozen-score audit rather
than a complete-panel replication or primary model validation.

GSE108215 is a separate, supplementary whole-blood MCI intervention-response
audit. Its official matrix contains 86 rows and 35 traceable pre/post subjects
(HEP=18, MAP=17). The GPL570 platform maps 10/12 exact frozen symbols; it has
no HC group and cannot validate the primary stage endpoint or establish an
intervention effect. The GEO narrative reports 51 pre-treatment rows, whereas
the processed matrix contains 50 pre-treatment and 36 post-treatment rows;
the audit uses only matrix-observed rows. GSE153881 remains blocked because its complete GEO sample
metadata and supplementary geneCOUNT filenames conflict on MCI/AD group
assignments.

The GSE1408 overlap audit records that GSE140831 contains all GSM records from
GSE140829 and GSE140830, while GSE140830 has no MCI samples. These series are
therefore excluded from independent MCI-versus-HC validation and are not
counted as additional cohorts.

## Reproduction order

```text
02_ANALYSIS\run_pairwise_external.py
02_ANALYSIS\run_study_heldout.py
02_ANALYSIS\run_leakage_attack.py
02_ANALYSIS\run_gse249477_stress_test.py
02_ANALYSIS\run_gse249477_rank_transport.py
02_ANALYSIS\run_tsp_relative_expression_transport.py
02_ANALYSIS\run_cross_study_stable_signal.py
02_ANALYSIS\run_domain_invariant_sensitivity.py
02_ANALYSIS\run_split_conformal.py
02_ANALYSIS\run_ood_distance_sensitivity.py
02_ANALYSIS\build_panel_provenance.py
02_ANALYSIS\audit_legacy_mr_method_boundary.py
02_ANALYSIS\audit_gse270454_join.py
02_ANALYSIS\run_gse270454_locked12_joined.py
02_ANALYSIS\run_gse108215_paired_audit.py
02_ANALYSIS\run_gse282742_progression_association.py
02_ANALYSIS\audit_gse1408_overlap.py
02_ANALYSIS\audit_gse63063_overlap.py
02_ANALYSIS\run_domain_separability_audit.py
02_ANALYSIS\audit_data_ledger.py
02_ANALYSIS\make_figure_pack.py
99_QA\run_full_qa.py
```

The primary and study-held-out outputs include an explicit sex-only baseline.
Calibration is reported as descriptive logistic calibration assessment after
predictions are locked; it is not used to recalibrate external scores. The external selective
risk curve is indexed only by development-derived confidence thresholds and is
restricted to the primary MCI-versus-HC endpoint. GSE140829 labels are parsed
from explicit `diagnosis:` metadata; title-derived subject tokens remain
candidate diagnostics because repeated and cross-diagnosis tokens do not
establish a validated subject map.

## Candidate discovery boundary

The saved official GEO query for `"mild cognitive impairment"` plus `blood`
returned 38 unique GSE records on 2026-08-13. Eighteen were absent from the
previous candidate registry; none met the prespecified human direct peripheral
whole-blood bulk-mRNA, independent MCI-versus-HC criterion. This is the scope of
one saved GEO query, not a statement that all public archives have been
exhausted. See `01_DATA_AUDIT/public_mci_blood_candidate_discovery_20260813.md`
and the saved NCBI query artifacts for the source metadata.

## Release status

The latest full QA report is `PASS` with zero failures. It checked 31 analysis
scripts, 7 main figures, 1 supplementary figure, and 46,300 cross-cohort effect
rows. The verified execution environment is Python 3.10.11 with the pinned
direct dependencies in `requirements-v5.txt`.

This directory is the curated V5 public release layer. The raw source matrices,
row-level prediction tables, subject-structure records and large uncompressed
cross-omics tables remain outside the repository; they are either retrievable
from the cited public accession or represented by aggregate audit tables.
Do not add raw participant-level data, sample identifiers, credentials, or
unverified external-validation claims to the public release.

The immutable repository tag for this release is `v5.0.0`. The repository DOI
will be added after the tag is archived.
