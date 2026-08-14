# Submission readiness assessment: V5

Date: 2026-08-13

## Working manuscript identity

The defensible paper is a leakage-aware reliability and transportability audit
of public blood-transcriptome MCI models. It is not a validated clinical
diagnostic biomarker paper and should not be submitted under that claim.

## Evidence-backed assessment

| Dimension | Current state | Editorial meaning |
|---|---|---|
| Scientific question | Clear: what survives leakage control and independent transport? | Stronger than a routine gene-list paper |
| Internal modelling | Fold-contained selection; nested AUC 0.795 | Supports a development result only |
| Same-family transfer | AUC 0.730-0.770 in bidirectional AddNeuroMed hold-out | Shows study-family transfer, not broad generalization |
| Independent transport | GSE140829 AUC 0.483; Brier 0.344; ECE 0.340 | Negative external result is the central limitation and also the audit finding |
| Domain-shift diagnosis | Domain-only AUC 1.000; permutation P=0.00020 | Strong methodological evidence of expression-domain displacement |
| Uncertainty | q95 OOD 0.859; conformal coverage 0.669 at nominal 0.90 | Supports a failure-boundary story |
| Frozen panel | 10/12 exact symbols in external cohorts | Blocks complete-panel external replication |
| Analysis unit | GSE140829 subject independence unresolved | Prevents subject-level clinical validation wording |
| Cohort independence | GSE63063 is exactly the 717-GSM aggregate of GSE63060/GSE63061 and is excluded | Prevents duplicate AddNeuroMed development data being counted as external evidence |
| Candidate compatibility | Three purified-neutrophil SubSeries and two BCG-linked CSF/PBMC single-cell studies are excluded | Official metadata show zero pairwise GSM overlap for the neutrophil SubSeries, but no subject identifiers and incompatible sex/APOE/IL18R1 selection; this preserves the whole-blood bulk MCI-versus-HC estimand and prevents false validation counts |
| Secondary cross-platform transport | GSE249477 rank/z AUC 0.576 (95% CI 0.387-0.776; n=41), Brier 0.488, ECE 0.488; no target calibration | Weak/uncertain assessment-only stress test, not positive replication or external validation |
| Relative-expression transport audit | Development-only TSP: internal AUC 0.767 (95% CI 0.694-0.836); GSE140829 AUC 0.548 (0.486-0.612); GSE249477 AUC 0.465 (0.300-0.629), with no target fitting or calibration | Relative expression does not stably restore transport; this is a negative method stress test, not external validation or clinical deployment |
| Cross-study stable-signal audit | Direction-concordant stable ridge: GSE63060/GSE63061 outer-fold AUC 0.840/0.779; 100-repeat full-pipeline permutation p=0.0099; GSE140829/GSE249477 AUC 0.541/0.524 | Stable disease-effect directions in a related development family do not restore external transport; the result improves the failure diagnosis, not the clinical validation level |
| Cross-cohort effect concordance | GSE63060/GSE63061 effect Spearman=0.673 and sign concordance=0.718; GSE249477 comparisons Spearman=-0.174 to -0.191 and sign concordance=0.426-0.444 | Adds a non-predictive, quantitative layer showing related-study agreement but cross-platform discordance; GSE140829 remains a sample-row audit and is excluded from pooled inference |
| Candidate discovery | Saved exact GEO MCI-plus-blood query returned 38 records; 18 newly screened and none eligible | Supports cohort-accounting boundary but is not global public-archive exhaustion |
| Reproducibility | Protocol, scripts, source data, hashes, QA and panel traceability present | Strong package-level reproducibility |
| Literature and frozen-panel provenance | Core reporting, external-validation, dataset-shift and conformal-method references are verified; GSE140829 is title-matched to a Crossref-registered bioRxiv preprint; the archived 95-to-83-to-12 development selection is source-hash-locked and replayed | The preprint link does not resolve subject independence; panel replay is development-only and does not reconstruct original raw preprocessing or establish external validation |

## Journal-level positioning

### CNS or four major medical journals

Current fit is low. The main barriers are the negative primary external AUC,
unresolved subject independence in the external audit, incomplete exact-symbol
mapping, absence of prospective or subject-level clinical utility validation, and
incomplete verified literature support. CNS-level figures improve readability
but cannot compensate for these evidence gaps. A CNS attempt would require a
substantially broader and independently replicated methodological contribution,
preferably with multiple subject-level sites and a demonstrated intervention or
clinical decision consequence.

### Upper-tier medical informatics, translational methods or reliability journals

Current fit is conditional but credible. The manuscript has a coherent negative
result plus methodological diagnosis story: internal discrimination and
same-family transfer are separated from independent transport, and severe domain
displacement is quantified without using disease labels for domain-model fitting.
The most plausible editorial audience is methodological, reproducibility,
biomedical informatics or biomarker-evaluation readers, rather than editors
seeking a new clinically validated MCI signature.

### Chinese Academy of Sciences Q2 or above

The package is a Q2 candidate in the right journal category, but Q2 status and
acceptance cannot be guaranteed. The probability is materially lower for
clinical biomarker journals that expect positive independent validation, and
higher for journals that explicitly publish methods, reproducibility,
transportability, dataset shift or trustworthy biomedical AI studies. Exact
impact factors and CAS quartiles are intentionally not listed because the local
workbook at `D:\投稿情况\影响因子` was unavailable during this run; these values
must be queried from the current workbook before selecting a journal.

## Hard blockers before submission

1. Complete retained supplementary-method source matches and document any
   remaining raw-preprocessing provenance limits. The frozen-panel archived
   development selection replay is now verified. GSE140829 is title-matched to the bioRxiv preprint
   `10.1101/2019.12.13.875112`; this is provenance only, not peer-reviewed
   publication or a validation upgrade. The core reporting and method records
   are listed in `source_verification_register.md`.
2. Resolve or explicitly retain the GSE140829 sample-to-subject ambiguity in
   the abstract, main text, tables and cover letter.
3. Keep all external frozen-panel results labelled partial sensitivity or
   boundary audit; never call them complete-panel replication.
4. Add the final public code repository, immutable release tag and archival DOI
   to Data Availability and Code Availability. The in-package direct-dependency
   pins and verified Python 3.10 execution record are available in
   `requirements-v5.txt` and `docs/reproducibility_environment_v1.md`; archive
   the resolved target-platform environment with the public release.
5. Decide whether the progression branch is supplementary only. It must not be
   used to imply predictive conversion validation without an independent
   subject-level progression cohort.
6. Recheck current target-journal author instructions and figure dimensions
   before submission.

7. Keep E-MTAB-6094 outside the primary endpoint: it is PBMC/mononuclear-cell
   AD-versus-control data with no MCI group and an official metadata count
   conflict. Its audit is a compatibility-gate result, not external validation.
8. Keep GSE63063 outside model fitting, held-out testing and independent-
   validation counts: it is the complete GSE63060/GSE63061 GSM aggregate.
9. Keep GSE243743/GSE243744/GSE243745 and GSE280720/GSE281606 outside the
   primary endpoint and all independent-validation counts; only prespecified,
   directional supplementary context is permissible after their unresolved
   subject-map or diagnosis-map limitations are retained. Zero GSM overlap among
   the neutrophil SubSeries does not establish participant independence.

## Main reviewer risks

- The domain AUC of 1.000 may reflect platform, batch, recruitment or other
  acquisition differences; the paper must not identify a biological cause.
- Sample-row bootstrap intervals may be optimistic if repeated external rows are
  present.
- GSE140829 candidate title prefixes cross HC/MCI/AD labels and are descriptive
  fields rather than a validated subject map. They must not be used for
  cluster-bootstrap intervals, deduplication or subject-level inference.
- A rank/z-score AUC of 0.605 in GSE140829 is a partial sensitivity result, not
  a solved domain-adaptation method. The GSE249477 secondary cross-platform
  audit is weak and uncertain (AUC 0.576; 95% CI 0.387-0.776) and has severe
  probability miscalibration (Brier/ECE 0.488); it cannot be presented as
  positive replication. The development-only relative-expression TSP transport
  audit did not solve this boundary (GSE140829 AUC 0.548; 95% CI 0.486-0.612;
  GSE249477 AUC 0.465; 95% CI 0.300-0.629).
- Computational MR, virtual knockout, LINCS, network, docking, MD, ADME and
  cell-level analyses remain hypothesis-generating and cannot rescue the primary
  transport result.
- The frozen 12-gene panel has a source-hash-locked, deterministic archived
  development-selection replay (95 candidates, 83 development-consistent,
  6 MCI-peak plus 6 MCI-trough) but lacks independent discovery and a complete
  reconstruction of the original raw-preprocessing lineage; it remains an audit
  comparator.
- Direction-concordant selection across two related AddNeuroMed series provides
  development-family consistency but is not a new independent study. Its locked
  GSE140829/GSE249477 AUCs (0.541/0.524) remain near chance; target labels were
  used only after prediction lock for descriptive assessment, not for mapping,
  selection, model fit, prediction calibration, thresholding or adaptation.
- Cross-cohort effect discordance is descriptive. It cannot identify whether
  platform, recruitment, cell composition, preprocessing or biology causes the
  observed pattern, and GSE140829 cannot enter subject-level or pooled inference.

## Current release state

The V5 release has seven PDF/SVG/600-dpi PNG figures, source data, a hash manifest,
panel-level traceability, a claim-evidence ledger, an English manuscript draft
and supplementary analysis boundaries. The development-only TSP transport audit
is integrated into Fig. 5, its source-data mirrors and Table S13 mapping. The
fresh 2026-08-13 release verification passed with 0 failures across 24
registered analysis scripts, 7 figures and 66 figure-manifest entries.

The manuscript now contains draft Data Availability and Code Availability
statements with verified official GEO links, a direct-dependency manifest, a
verified Python 3.10 execution record and a runtime smoke test. These local
artifacts make the package ready to archive, but they do not create a public
repository, immutable tag, archive license or DOI.

These are package verification results, not evidence of journal acceptance.

## New candidate-screening result

The official GEO metadata audit of GSE140829, GSE140830 and GSE140831 found
that GSE140831 contains all 587 GSM records from GSE140829 and all 542 GSM
records from GSE140830. GSE140830 has no MCI samples. These records are thus
not additional independent validation cohorts and are excluded from the
primary MCI-versus-HC endpoint. This strengthens the audit's independence
accounting but does not increase the number of external cohorts.

The official BioStudies audit of E-MTAB-6094 found 35 factor-table samples
(normal=13, AD CDR1=15, AD CDR3=7), no MCI samples, PBMC/mononuclear-cell input,
and a conflict with the study description's 25 AD plus 15 controls. The series
was excluded from primary validation and no grouped expression analysis was
performed.

The cached official GEO summary also identifies GSE63063 as the complete
717-GSM union of GSE63060 (329 GSMs) and GSE63061 (388 GSMs), with no additional
GSE63063 records. It is therefore an accession-level aggregate audit, not a
third AddNeuroMed cohort, and is excluded from model inputs and independent
validation. Three same-family purified-neutrophil SubSeries
(GSE243743/GSE243744/GSE243745) and two BCG-linked CSF/PBMC single-cell studies
(GSE280720/GSE281606) are registered as incompatible with the whole-blood bulk
MCI-versus-HC estimand. Automated contract tests require all five to remain
excluded from the domain registry and data ledger.

The saved official GEO query for `"mild cognitive impairment"` plus `blood`
returned 38 unique GSE records on 2026-08-13. Eighteen records were not in the
prior registry; none met the prespecified criteria for a human, direct
peripheral whole-blood bulk-mRNA, independent MCI-versus-HC comparator. This is
an explicit GEO-query result, not a claim that all public repositories have been
exhausted. The discovery note and raw metadata artifacts are retained in the
release package.
