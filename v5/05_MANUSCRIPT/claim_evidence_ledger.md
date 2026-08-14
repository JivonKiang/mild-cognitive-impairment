# Claim-evidence ledger for the MCI public-data manuscript

Status: working manuscript boundary, V5 generated 2026-08-12.

## Evidence tiers

- **Supported**: public source data, inspectable code, and result artifact are all present.
- **Exploratory**: result is available but limited by mapping, sample size, platform, or lack of independent validation.
- **Blocked**: the proposed claim cannot be made from the current evidence.

## Legacy MR exclusion audit

**Blocked claim:** inherited cis-eQTL-to-AD MR or colocalization output
establishes causal MCI biology, a validated target, therapeutic efficacy or a
clinical biomarker.

**Evidence:** `02_ANALYSIS/audit_legacy_mr_method_boundary.py` and
`03_RESULTS/legacy_mr_method_audit.json`. The inherited outputs contain 85 MR
rows, including 8 FDR-positive rows, but no colocalization result with
`PP.H4 >= 0.5`. More importantly, the code-and-provenance audit fails V5 gates
for recorded exposure effect scale, LD clumping, gene-locus-specific
colocalization input, MCI-stage phenotype compatibility and an inspectable
rerun environment.

**Allowed wording:** “Inherited MR/colocalization output was excluded from
causal interpretation after a code-and-provenance audit.”

**Forbidden wording:** “MR proved causality,” “colocalization validated the
target,” “causal MCI gene,” or “therapeutic target.”

| Proposed claim | Tier | Evidence | Allowed manuscript wording | Forbidden wording |
|---|---|---|---|---|
| Fold-level feature selection prevents a known leakage route | Supported | `run_pairwise_external.py`; `leakage_attack_summary.csv`; repeated 5-fold diagnostic | “The prespecified fold-contained pipeline reduced a label-aware selection route and separated internal optimism from transportability.” | “The pipeline eliminates all possible bias.” |
| Transcriptome model discriminates MCI from HC in the development family | Supported, internal only | `metrics.csv`, nested internal row; n=208 | “Within the AddNeuroMed development family, nested internal AUC was 0.795.” | “The model is clinically accurate.” |
| Model externally validates in GSE140829 | Blocked as a positive validation claim | `metrics.csv`; AUC=0.483, CI 0.420–0.542; OOD=0.859 | “Performance did not transport to the independent GSE140829 cohort.” | “The signature was externally validated.” |
| GSE140829 title prefixes support cluster bootstrap, deduplication or subject-level inference | Blocked | `gse140829_subject_structure_audit.json`; candidate prefixes cross HC/MCI/AD labels and are descriptive title fields only | “GSE140829 remains a sample-row transport audit because a validated sample-to-subject map is unavailable.” | “Candidate title prefixes define subjects,” “cluster-bootstrap confidence intervals,” or “subject-level external validation.” |
| Same-platform study-held-out transfer is possible | Exploratory | `study_heldout_metrics.csv`; bidirectional AUC 0.730–0.770 | “Transfer between the two AddNeuroMed series was higher than transfer to GSE140829, indicating study-family dependence.” | “The model generalizes across populations.” |
| Frozen 12-gene panel is a stable biomarker | Blocked | Exact mapping blocked in external cohorts; partial AUC 0.580 and 0.410 | “The locked score was retained as an audit comparator; external exact-symbol mapping was incomplete.” | “The 12-gene panel is a validated biomarker.” |
| GSE249477 independently replicates the MCI signal | Blocked | `gse249477_stress_test.json`; n=20 MCI/21 HC; partial 10/12; AUC=0.410 | “GSE249477 was used as a cross-platform stress test and did not provide positive replication.” | “GSE249477 confirmed the signature.” |
| Progression signature predicts conversion | Blocked | GSE282742 progression branch; n=52 unique baseline subjects; no external progression cohort; source artifact must be linked before manuscript release | “Progression analyses were exploratory and underpowered for predictive validation.” | “The signature predicts future AD conversion.” |
| GSE282742 Reactome enrichment establishes progression mechanism | Exploratory only | `gse282742_progression_pathway_summary.csv`; 20,897 dual-source-mapped symbols; 1,396 tested pathways; 47 BH-FDR <0.05 | "In one small public representative-record cohort, adjusted P-MCI/S-MCI rank differences showed exploratory Reactome pathway associations." | "The pathways establish cell-type changes, MCI-to-AD mechanism, causal biology, or a conversion biomarker." |
| MR, virtual knockout, LINCS, network, docking or MD establish mechanism | Blocked | Supplementary computational outputs only | “These analyses nominate hypotheses for future experimental testing.” | “These analyses demonstrate mechanism or therapeutic efficacy.” |
| The public-data workflow is clinically deployable | Blocked | no prospective cohort, no clinical utility validation, strong OOD | “The findings define a reliability boundary for public blood-transcriptome biomarker claims.” | “The model is ready for clinical use.” |
| Label-free sample-wise rank/z-score resolves domain shift | Exploratory | `domain_invariant_sensitivity.csv`; external AUC 0.605, calibration remains limited | “A label-free sample-wise transform partly improved ranking discrimination but did not establish calibrated transportability.” | “Domain invariance solved external validation.” |
| Target-unlabeled CORAL alignment restores transportability | Blocked as a transportability or deployment claim | `run_transductive_coral_audit.py`; `transductive_coral_metrics.csv`; `transductive_coral_audit.json`; source-fitted preprocessing/PCA; target-unlabeled latent CORAL; GSE140829 AUC 0.544 (95% CI 0.484-0.597), Brier 0.276, ECE 0.198; GSE249477 AUC 0.421 (95% CI 0.255-0.610), Brier 0.294, ECE 0.187 | “A target-unlabeled, batch-level CORAL stress test did not establish calibrated transportability in either external assessment.” | “CORAL externally validated the model, solved domain shift, supports single-sample deployment, or established a generalizable MCI biomarker.” |
| ComBat is a valid label-free alternative in the present targets | Blocked as a method-eligibility claim | ComBat source: DOI 10.1093/biostatistics/kxj037. Target MCI/HC diagnoses are intentionally unavailable before assessment, so a disease covariate cannot be protected during batch correction. | “ComBat was not applied because a target-label-free adjustment could not preserve the primary disease covariate under this design.” | “ComBat was tested and failed”, or “ComBat can be safely applied without target disease covariates.” |
| Development-only relative-expression TSP restores cross-platform MCI discrimination | Blocked as a transportability claim | `tsp_transport_metrics.csv`; `tsp_transport_selected_pairs.csv`; `tsp_transport_audit.json`; internal AUC 0.767, GSE140829 AUC 0.548 (95% CI 0.486-0.612), GSE249477 AUC 0.465 (95% CI 0.300-0.629) | “A development-only relative-expression TSP transport audit retained internal discrimination but did not stably restore transport across the two external assessments.” | “TSP externally validated an MCI biomarker, solved domain shift, or is clinically deployable.” |
| Direction-concordant cross-study stable signal restores external transport | Blocked as a transportability claim | `cross_study_stable_metrics.csv`; `cross_study_stable_selected_genes.csv`; `cross_study_stable_permutation_null.csv`; AddNeuroMed outer-fold AUC 0.840/0.779, permutation p=0.0099, GSE140829 AUC 0.541 (95% CI 0.479-0.598), GSE249477 AUC 0.524 (95% CI 0.500-0.575) | “The direction-concordant model demonstrated development-family consistency but did not restore locked external transport.” | “The stable signal is independently validated, clinically calibrated, or solves domain shift.” |
| OOD failure is robust to distance definition | Exploratory | `ood_distance_sensitivity.csv`; mean-absolute standardized OOD 0.574-0.963 across q99-q90; Ledoit-Wolf and robust PCA OOD 1.00 across tested thresholds | “Severe external domain displacement was observed across training-only distance definitions, although the absolute OOD rate depends on the detector and threshold.” | “The OOD detector provides a validated probability of domain membership.” |
| Split-conformal prediction provides valid external coverage | Blocked for target-domain coverage | `split_conformal_summary.csv`; mean external coverage 0.669 at nominal 0.90 | “Development-calibrated conformal sets were evaluated as a stress test; nominal coverage did not transfer to the external domain.” | “Conformal prediction guarantees clinical coverage across domains.” |
| Expression domains are separable after phenotype metadata are discarded from model fitting | Supported as domain-shift evidence | `domain_separability_summary.csv`; domain-only AUC=1.000, repeated 5-fold, permutation p=0.00020 | “Development and external expression domains were perfectly separable in this audit after phenotype metadata were discarded from model fitting, demonstrating severe domain displacement.” | “The domain classifier proves the biological cause of disease-model failure.” |
| Gene-level MCI-HC effects are stable across all public cohorts | Blocked as a universal-stability claim | `cross_cohort_effects.csv`; `cross_cohort_effect_pairwise_summary.csv`; `cross_cohort_effect_heterogeneity.csv`; AddNeuroMed Spearman=0.673 and sign concordance=0.718, whereas GSE249477 comparisons were negative (Spearman -0.174 to -0.191) | “The effect audit found substantial within-AddNeuroMed directional agreement but materially discordant cross-platform effects in GSE249477; GSE140829 is a sample-row audit with unresolved subject independence.” | “All cohorts replicated a common blood MCI transcriptional effect,” or “observed discordance proves a biological mechanism.” |

| GSE108215 supports MCI diagnostic replication or intervention efficacy | Blocked | `gse108215_paired_audit.json`; all participants have MCI; 35 paired subjects; partial 10/12 score; no cognitive endpoint analysis | A small public MCI intervention cohort provided exploratory paired molecular-response context for the frozen partial score. | The signature replicated MCI diagnosis, or MAP improved cognition / was clinically effective. |
| GSE153881 supports plasma-EV MCI/AD stress-test results | Blocked pending label reconciliation | Official series matrix and complete RAW filenames invert MCI and AD accession assignments | GSE153881 was registered as a candidate but no grouped analysis was performed because its public label sources conflict. | GSE153881 replicated or refuted the MCI signature. |
| GSE140830 or GSE140831 add independent blood validation | Blocked as independent replication | `gse1408_overlap_audit.json`; GSE140830 has no MCI samples; GSE140831 shares all 587 GSE140829 GSM records and all 542 GSE140830 GSM records | These accessions were excluded from independent validation and are not counted as additional cohorts. | GSE140830/831 independently replicated the MCI model. |
| E-MTAB-6094 adds an independent MCI-HC blood validation cohort | Blocked | `emtab6094_compatibility_audit.json`; official BioStudies factor table has 13 normal, 15 AD CDR1, 7 AD CDR3 and 0 MCI samples; specimen is mononuclear cells/PBMC; description and factor-table totals conflict | E-MTAB-6094 was excluded from the primary endpoint and retained only as an orthogonal PBMC/AD-HC compatibility-gate example pending SDRF reconciliation. | E-MTAB-6094 independently validated the whole-blood MCI model. |
| GSE63063 adds a new AddNeuroMed external or held-out validation cohort | Blocked | `gse63063_overlap_audit.json`; all 329 GSE63060 and all 388 GSE63061 official GSM records occur in GSE63063; aggregate union=717 and zero GSE140829 GSM overlap | GSE63063 is an accession-level aggregate audit and is excluded from model fitting, held-out testing and independent-validation counts. | GSE63063 independently replicated, enlarged, or externally validated the AddNeuroMed model. |
| GSE243743/GSE243744/GSE243745 validate the whole-blood MCI-HC model | Blocked | Candidate registry; official GEO records; [PMID 38961225](https://pubmed.ncbi.nlm.nih.gov/38961225/) | These same-family, purified-neutrophil records may be discussed only as conditional cell-type-specific context after a subject-map audit. | The neutrophil SubSeries are independent whole-blood validation, a pooled external cohort, or causal proof. |
| GSE208623 validates the blood mRNA model or establishes methylation mechanism | Blocked | Source-locked GEO beta matrix and EPIC annotation; sex-adjusted MCI-HC association in 20 MCI/20 HC rows; 20 AD rows excluded; 10,422 gene overlap; age and cell composition unadjusted | GSE208623 provides an independent-source, cohort-level orthogonal cross-omics association audit. It does not upgrade the negative mRNA transport result. | GSE208623 validated/tuned/calibrated the mRNA model, replicated the signature, established methylation-expression mediation, proved cell composition, or demonstrated causal biology. |
| GSE208623 has a descriptive cross-omics association with the locked transcriptome reference | Exploratory | `gse208623_cpg_associations.csv`; `gse208623_gene_concordance.csv`; `gse208623_cross_omics_source_audit.json`; 738,242 CpGs; 6,869 BH FDR <0.05; Spearman rho=0.280; same-sign fraction=0.560; fixed 10,000-permutation P=0.00010 | The source-locked, sex-adjusted audit found descriptive gene-level concordance between GSE208623 methylation association aggregates and the frozen AddNeuroMed transcriptome effects. | Independent gene-level differential-methylation testing, external validation, biological mechanism, mediation, causal inference, or clinically useful prediction. |
| GSE280720/GSE281606 validate a blood MCI-HC model or BCG mechanism | Blocked | Candidate registry; official GEO records; [NCT05004688](https://clinicaltrials.gov/study/NCT05004688); [PMID 42393277](https://pubmed.ncbi.nlm.nih.gov/42393277/) | These longitudinal CSF/PBMC single-cell records are intervention-linked immune-cell context only, with pairing and exposure retained. | They provide MCI-HC discrimination, blood-bulk replication, disease mechanism, or vaccination efficacy evidence. |

## GSE249477 rank/z transport boundary

The separately locked rank/z sensitivity in GSE249477 is an exploratory, negative or
uncertain secondary independent cross-platform transport audit. Its results are
recorded in `gse249477_rank_transport_metrics.csv` and
`gse249477_rank_transport_audit.json`: 41 sample rows (21 HC and 20 MCI), AUC
0.576 (bootstrap 95% CI 0.387-0.776), Brier 0.488 and ECE 0.488. The permitted
wording is that ranking was weak and uncertain and probability calibration was
poor. It is forbidden to call this external validation, positive replication,
a solution to domain shift, or evidence of a generalizable MCI blood biomarker.

The saved 2026-08-13 GEO query and candidate screen are recorded in
`public_mci_blood_candidate_discovery_20260813.md`. This 38-record query did
not identify a previously unregistered independent direct whole-blood bulk-mRNA
MCI-HC cohort meeting the prespecified criterion. It is forbidden to describe
this explicit GEO query as exhaustion of all public archives.

## Development-only relative-expression TSP transport boundary

The separately locked TSP extension is a development-only relative-expression
transport audit, motivated by the pairwise mRNA comparison method family
([PMID 16646797](https://pubmed.ncbi.nlm.nih.gov/16646797/) and [PMID
23496976](https://pubmed.ncbi.nlm.nih.gov/23496976/)). The methods literature
only supports testing this classifier family; it does not support MCI clinical
utility or an expectation of positive transport. Within the AddNeuroMed
development rows, five-fold out-of-fold AUC was 0.767 (95% CI 0.694-0.836).
After the final ten disjoint pairs were selected exclusively in AddNeuroMed,
GSE140829 scored 0.548 (95% CI 0.486-0.612; n=383) and GSE249477 scored 0.465
(95% CI 0.300-0.629; n=41). All 587 GSE140829 rows and all 62 GSE249477 rows
were predicted before endpoint assessment; GSE249477 diagnosis labels were not
parsed until after target scoring, and no target calibration was fitted. The
audit does not establish external validation, positive replication, target
calibration, clinical utility, or a solution to domain shift.

## Cross-study stable-signal transport boundary

The cross-study stable-signal extension is a public-data transport stress test,
not an independent-discovery or clinical-biomarker study. Before any target
diagnosis metadata were read, the candidate universe was fixed to 11,575 exact
symbols shared by GSE63060, GSE63061, GSE140829 and GSE249477. In each
synchronized AddNeuroMed outer fold, MCI-HC effect directions were calculated
only from the two study-specific training partitions. Non-zero, same-direction
genes were deterministically ranked by their minimum absolute study effect and a
fixed 100-gene ridge was fit on the combined training rows. The resulting AUCs
were 0.840 in GSE63060 and 0.779 in GSE63061; the pooled development OOF AUC
was 0.819 with a 100-repeat full-pipeline label-permutation p=0.0099.

All 587 GSE140829 and all 62 GSE249477 expression rows were scored before their
disease labels were used for assessment. Locked MCI-HC AUC was 0.541 (95% CI
0.479-0.598; n=383) in GSE140829 and 0.524 (95% CI 0.500-0.575; n=41) in
GSE249477. Target labels were not used for mapping, selection, preprocessing,
fitting, prediction calibration, thresholding or adaptation; after predictions
were locked, they were used only for descriptive performance and calibration
assessment. This result supports development-family directional consistency but
does not establish independent validation, calibrated transport, clinical
utility or a solution to domain shift.

## Cross-cohort effect heterogeneity boundary

The non-predictive cross-cohort effect audit fixes the 11,575-gene exact-symbol
universe before external diagnosis labels are accessed and computes independent
MCI-minus-HC Hedges' g within each cohort. GSE63060 and GSE63061 show a
Spearman effect correlation of 0.673 and 71.8% sign concordance across 11,575
estimable genes. The AddNeuroMed unweighted descriptive reference has a
Spearman correlation of 0.486 and 66.2% sign concordance with GSE140829.
By contrast, GSE249477 comparisons use 10,923 estimable genes because 652
shared genes are non-estimable there; their effect correlations with the other
three cohorts range from -0.174 to -0.191 and sign concordance from 42.6% to
44.4%.

These are descriptive effect-concordance results, not pooled estimates or
cross-cohort heterogeneity tests. The two-study Q/I-squared output is restricted
to GSE63060/GSE63061. GSE140829 remains excluded from pooled inference because
subject independence is unresolved. The audit does not prove assay, cell-type
or biological causality; it quantifies an additional, label-safe reliability
boundary for general MCI blood-transcriptome claims.

## Mandatory data corrections

- GSE270454 count-column to GEO metadata join: 45/45 unique whole-blood RNA
  samples after mapping descriptive cohort labels and terminal sample ordinals.
- GSE270454 frozen-score boundary: 10/12 exact symbols; MCI-HC AUC=0.404
  (bootstrap 95% CI 0.204-0.612); age/sex-adjusted estimate=-0.147,
  P=0.373; n=45 rows.

- GSE249477 official series titles: HC=21, MCI=20, AD=21; total n=62.
- GSE249477 rank/z transport audit: 62/62 title-to-TPM identifiers matched and
  all 62 target rows were transformed before target diagnosis parsing. The
  secondary HC/MCI assessment contains 41 sample rows (21 HC, 20 MCI); the 21 AD
  rows are excluded only after transformation. GSE249477 labels are false for
  target transformation, fitting, feature selection and thresholding/calibration
  and true only for the locked final assessment. No target calibration is fit.
- GSE282742: 116 sample rows, 86 unique subjects, 52 metadata-order representative P-MCI/S-MCI records (22 P-MCI, 30 S-MCI), 51 age/sex-complete records (22 P-MCI, 29 S-MCI), and 11 representative P-MCI-to-AD paired subjects. Public metadata provides no verified visit dates.
- The 52-subject progression analysis and the 116-row cohort audit must never be reported as the same `n`.
- GSE108215: 86 expression rows include 35 traceable paired MCI subjects (HEP=18;
  MAP=17). The exact frozen panel maps 10/12 GPL570 symbols; `PHF15` and
  `LAT1-3TM` are absent. It is supplementary paired molecular-response context,
  not stage validation or evidence of intervention efficacy.
- GSE108215 cohort accounting: the official series narrative reports 51
  pre-treatment and 36 post-treatment samples, while the processed expression
  matrix contains 50 pre-treatment and 36 post-treatment columns. All paired
  summaries use the observed 86 matrix rows only.
- GSE153881: the public GEO sample metadata and 15 complete supplementary
  `geneCOUNT` filenames reverse MCI/AD accession assignments. The series is
  blocked pending independent authoritative label reconciliation; no groupwise
  result is reported.
- GSE1408 overlap audit: GSE140831 contains all GSM records in GSE140829 and
  GSE140830; GSE140830 contains no MCI samples. Neither accession is counted as
  independent MCI-versus-HC validation.
- E-MTAB-6094: official BioStudies metadata describe blood-derived mononuclear
  cells/PBMC, not whole blood. The factor table contains 35 samples (13 normal,
  15 AD CDR1 and 7 AD CDR3) and no MCI, whereas the study description states 25
  AD and 15 controls. It is excluded from the primary endpoint pending SDRF
  reconciliation; no grouped expression analysis was performed.

## Updated primary external audit boundary

GSE140829 metadata provide explicit `diagnosis:` labels: HC=249, MCI=134,
AD=204. Title-derived candidate subject tokens repeat across 257 sample rows
and cross diagnosis labels; they are not a validated subject map. The primary
external results therefore remain sample-row transport audits, not
subject-level validation claims. The data ledger uses the same role and unit
labels: external transport assessment and sample-row audit.

The primary baseline comparison is: age-only AUC 0.518, sex-only AUC 0.452,
age-plus-sex AUC 0.498, and transcriptome AUC 0.483 in GSE140829. The
transcriptome model did not add external discrimination beyond the clinical
baselines.

GSE270454 has a verified count-column to GEO metadata join for 45/45
whole-blood RNA samples. Because the frozen panel maps only 10/12 exact
symbols and the MCI subgroup contains 10 rows, it is retained as an
independent small-cohort boundary audit. The partial score produced a
row-level MCI-HC AUC of 0.404 (bootstrap 95% CI 0.204-0.612); the descriptive
age/sex-adjusted MCI-HC estimate was -0.147 (P=0.373). These results are not
complete-panel replication or clinical validation.

Calibration intercept and slope in the metrics table are descriptive
evaluation-set assessments fit after predictions were locked. They are not
external recalibration parameters and were not used to alter external scores.

## Supplementary-only analyses

MR, virtual knockout, LINCS perturbation, network pharmacology, docking, molecular dynamics, ADME and cell-level simulations remain supplementary hypothesis-generation. They must be explicitly marked as computational and non-causal, and each panel must link to its own source data and code.

The V5 provenance ledger records the frozen panel's historical source table,
direction sign, fixed audit weight, exact-symbol mapping status and missing
symbol boundary. The archived 95-candidate to 83 development-consistent to
6-plus-6 locked-panel selection is deterministically replayed and source-hash
locked in `archived_locked12_selection_audit.*`. This is a development-only
replay across the related GSE63060/GSE63061 system, not independent discovery,
external validation, clinical biomarker validation, or a reconstruction of the
original raw-preprocessing lineage; the panel remains an audit comparator.
