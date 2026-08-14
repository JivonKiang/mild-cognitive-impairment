# Source Verification Register

Status date: 2026-08-14. Official database records were checked online where
marked `VERIFIED_ONLINE`. Original cohort-paper matching remains separate from
accession-level metadata verification.

| Source family | Record | Status | URL or identifier | Boundary |
|---|---|---|---|---|
| Reporting | TRIPOD statement | VERIFIED_ONLINE | PMID 25560714; DOI 10.7326/M14-0697; https://pubmed.ncbi.nlm.nih.gov/25560714/ | Reporting standard, not validation evidence |
| Reporting | TRIPOD+AI statement | VERIFIED_ONLINE | PMID 38626948; DOI 10.1136/bmj-2023-078378; https://pubmed.ncbi.nlm.nih.gov/38626948/ | Reporting standard, not validation evidence |
| Bias appraisal | PROBAST | VERIFIED_ONLINE | PMID 30596875; DOI 10.7326/M18-1376; https://pubmed.ncbi.nlm.nih.gov/30596875/ | Risk-of-bias/applicability tool |
| Observational reporting | STROBE | VERIFIED_ONLINE | PMID 18064739; DOI 10.1016/S0140-6736(07)61602-X; https://pubmed.ncbi.nlm.nih.gov/18064739/ | Reporting standard |
| Prediction validation | Internal, internal-external and external validation | VERIFIED_ONLINE | DOI 10.1016/j.jclinepi.2015.04.005; https://doi.org/10.1016/j.jclinepi.2015.04.005 | Supports separating validation roles |
| Dataset shift | A unifying view on dataset shift in classification | VERIFIED_ONLINE | DOI 10.1016/j.patcog.2011.06.019; https://doi.org/10.1016/j.patcog.2011.06.019 | Conceptual/method source; does not identify a biological cause |
| Conformal inference | Distribution-Free Predictive Inference for Regression | VERIFIED_ONLINE | DOI 10.1080/01621459.2017.1307116; https://doi.org/10.1080/01621459.2017.1307116 | Exchangeability assumptions do not imply target-domain validity under shift |
| Target-unlabeled domain adaptation | Return of Frustratingly Easy Domain Adaptation (CORAL) | VERIFIED_ONLINE | DOI 10.1609/aaai.v30i1.10306; https://ojs.aaai.org/index.php/AAAI/article/view/10306; accessed 2026-08-14 | Supports an unlabeled-target second-order alignment stress-test method only; it does not establish biomedical transport, clinical utility, or single-sample deployment. |
| Batch adjustment | Adjusting batch effects in microarray expression data using empirical Bayes methods (ComBat) | VERIFIED_ONLINE | DOI 10.1093/biostatistics/kxj037; https://doi.org/10.1093/biostatistics/kxj037; accessed 2026-08-14 | The method is not applied in the target-label-free audit because preserving the primary disease covariate would require unavailable target diagnostic covariates. |
| Relative-expression classifier method | Pairwise mRNA comparison classifier | VERIFIED_ONLINE | PMID 16646797; https://pubmed.ncbi.nlm.nih.gov/16646797/; accessed 2026-08-13 | Supports the TSP/pairwise-relative-expression method family only; not MCI clinical utility, positive transport, or deployment evidence |
| Relative-expression classifier method | AUREA relative-expression methods | VERIFIED_ONLINE | PMID 23496976; https://pubmed.ncbi.nlm.nih.gov/23496976/; accessed 2026-08-13 | Supports the relative-expression classifier method family only; not MCI clinical utility, positive transport, or deployment evidence |
| Pathway resource | Reactome current pathway GMT | VERIFIED_ONLINE | https://reactome.org/download/current/ReactomePathways.gmt.zip; accessed 2026-08-13; archive SHA-256 `8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4` | Supports the source-locked supplementary pathway-membership resource only; not MCI mechanism or clinical utility. |
| Gene-identifier cross-check | MyGene.info human Ensembl-gene symbol query | VERIFIED_ONLINE | https://mygene.info/v3/query; accessed 2026-08-13; local query snapshot in `01_DATA_AUDIT/reference/gse282742_mygene_ensembl_symbol_snapshot_20260813.json` | Supports the dual-source mapping gate for the GSE282742 supplementary analysis only; mismatched or unmapped IDs are excluded. |
| Official dataset record | GSE63060 | ACCESSION_AND_LINKED_PAPER_VERIFIED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060; PMID 26343147; DOI 10.1186/s13059-015-0750-x | AddNeuroMed development-family batch 1; not independent external validation |
| Official dataset record | GSE63061 | ACCESSION_AND_LINKED_PAPER_VERIFIED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061; PMID 26343147; DOI 10.1186/s13059-015-0750-x | AddNeuroMed development-family batch 2; not independent external validation |
| Official dataset record | GSE63063 | CACHED_OFFICIAL_METADATA_AUDITED; LIVE_RECHECK_PENDING_NETWORK_ACCESS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63063 | Complete GSM aggregate of GSE63060/GSE63061 (717 official GSM records); excluded from model and independent-validation roles |
| Official dataset record and title-matched preprint | GSE140829 | ACCESSION_AND_TITLE_MATCHED_PREPRINT_VERIFIED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140829; [10.1101/2019.12.13.875112](https://doi.org/10.1101/2019.12.13.875112); Crossref record https://api.crossref.org/works/10.1101/2019.12.13.875112; accessed 2026-08-14 | Official GEO title exactly matches the bioRxiv preprint title. This is a preprint linkage, not a peer-reviewed-paper claim; whole-blood array (GPL15988) remains a sample-row transport audit with subject independence unresolved. |
| Official dataset record | GSE249477 | ACCESSION_AND_LINKED_PAPER_VERIFIED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249477; PMID 40034360; DOI 10.1177/25424823241307878 | Small independent RNA-seq cross-platform stress test; not positive replication or clinical validation |
| Official dataset record | GSE282742 | LINKED_PAPER_VERIFIED_ACCESSION_HEADER_NOT_LOCAL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282742; PMID 39754267; DOI 10.1186/s13195-024-01651-0 | Separate progression-association branch; not independent MCI-vs-HC predictive validation or causal evidence |
| Official dataset record | GSE108215 | ACCESSION_AND_LINKED_PAPER_VERIFIED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108215; PMID 29900257; DOI 10.1016/j.dib.2018.03.086 | Supplementary paired intervention audit only; no HC comparator |
| Official dataset record | GSE270454 | ACCESSION_AND_LINKED_PAPER_VERIFIED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270454; PMID 40779089 | Small frozen-score boundary audit |
| Official BioStudies record | E-MTAB-6094 | FIELD_LEVEL_RECORD_VERIFIED; RAW_RESPONSE_NOT_PERSISTED | https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-6094 | PBMC/AD-HC, no MCI, count conflict; excluded from primary endpoint |
| Official GEO records | GSE243743/GSE243744/GSE243745/GSE243750 | ACCESSION_METADATA_REVIEWED; SUBJECT_MAP_PENDING | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243743; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243744; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243745; PMID 38961225 | Same-study-family purified-neutrophil context; not whole blood and not independent validation |
| Official GEO record | GSE208623 | ACCESSION_METADATA_AND_PROCESSED_BETA_SOURCE_LOCKED | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208623; PMID 37388929; processed beta endpoint https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208623/suppl/GSE208623_Norm_beta.csv.gz; beta SHA-256 `bf6c3f81bb676008a6c02630858b95a481ff5842200ee75f260146907a3500b3`; accessed 2026-08-14 | Chinese peripheral-blood-leukocyte EPIC methylation cohort, 20 HC/20 MCI/20 AD. The source-locked beta file supports the separate sex-adjusted cross-omics association audit only; it is not mRNA validation, model fitting/tuning/calibration, mediation, or causal evidence. |
| Literature-reported archive candidate | HRA019152 | LITERATURE_REPORTED_ARCHIVE; DIRECT_RECORD_METADATA_PENDING | PMID 42529199; PMCID PMC13416258; https://pmc.ncbi.nlm.nih.gov/articles/PMC13416258/; reported GSA-Human accession HRA019152; accessed 2026-08-14 | Article reports a 48-participant Chinese PBMC bulk-RNA-seq study with SCD/MCI/dementia each split by A-beta PET status. PBMC input, no healthy-control group and eight MCI rows per PET arm exclude it from the whole-blood MCI-HC estimand; its archive access state is not inferred from the article. |
| Official GEO accession mismatch audit | GSE165090 | OFFICIAL_ACCESSION_CONTRADICTS_SECONDARY_DATA_AVAILABILITY_CITATION | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165090; cited as AddNeuroMed ANM1 in PMID 36031604 / PMCID PMC9422129; accessed 2026-08-14 | Official GEO identifies GSE165090 as a trophoblast SERPINA1 RNA-seq study, not a blood MCI cohort. This documents an accession-verification failure mode; it is not a candidate cohort or biological result. |
| Official GEO/registry records | GSE280720/GSE281606 | ACCESSION_METADATA_REVIEWED; DIAGNOSIS_MAP_AND_OVERLAP_PENDING | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280720; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281606; https://clinicaltrials.gov/study/NCT05004688; PMID 42393277 | BCG-linked CSF/PBMC longitudinal single-cell context; excluded from primary endpoint |
| Official GEO discovery query | `"mild cognitive impairment" AND blood AND GSE` | SAVED_OFFICIAL_EUTILITIES_METADATA | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds | 38 unique GSE records returned on 2026-08-13; 18 previously unregistered records screened; explicit GEO query scope only, not global archive exhaustion |
| Frozen panel | Archived 95-candidate to 83 development-consistent to 12-gene selection chain | ARCHIVED_DEVELOPMENT_SELECTION_REPLAY_VERIFIED | `02_ANALYSIS/audit_archived_locked12_selection.py`; `03_RESULTS/archived_locked12_selection_audit.*`; archived scripts and tables SHA-256 locked | Deterministic GSE63060/GSE63061 development-only replay; not independent discovery, external validation, clinical biomarker validation, or original raw-preprocessing reconstruction |
| Legacy cis-eQTL-to-AD MR and colocalization | inherited R code and aggregate output tables | EXCLUDED_FROM_CAUSAL_INTERPRETATION_LEGACY_METHOD_AUDIT | `02_ANALYSIS/audit_legacy_mr_method_boundary.py`; `03_RESULTS/legacy_mr_method_audit.json`; method-source verification pending network recovery | The code-and-provenance audit identifies V5 gate failures and excludes inherited results from causal, MCI-stage, target-validation and biomarker claims. This does not judge correctly implemented MR or colocalization methods. |
| Public genetic-triangulation feasibility | GWAS Catalog metadata and eQTL Catalogue availability check | VERIFIED_ONLINE_NO_ANALYSIS_EXECUTED | `01_DATA_AUDIT/genetic_triangulation_feasibility_20260814.md`; GWAS Catalog metadata release 2026-08-02; MCI EFO query returned zero records; AD study `GCST90444373` has public full summary statistics; eQTL Catalogue `v3/studies` returned HTTP 500 | AD may be a future cross-stage context only. No MCI outcome or source-verified blood eQTL exposure was admitted, so no MR or colocalization was run. |
| Other supplementary computational analyses | virtual knockout, LINCS, network, docking, MD, ADME, cell simulation | PENDING_SOURCE_MATCH | `05_MANUSCRIPT/supplementary_analysis_index.md` | Hypothesis generation only |

The GSE208623 source-locked cross-omics audit is reproducible from
`02_ANALYSIS/run_gse208623_cross_omics_audit.py`; the privacy-safe aggregate
metadata audit is `02_ANALYSIS/audit_gse208623_metadata.py`. The source-locked
analysis is explicitly separate from the mRNA transport workflow.

The HRA019152 and GSE165090 candidate-screening extension is recorded in
`01_DATA_AUDIT/public_mci_blood_candidate_discovery_20260814.md`. It adds no
model input, fitted result, or validation cohort.

The E-MTAB-6094 field-level audit is reproducible from
`02_ANALYSIS/audit_emtab6094_compatibility.py`. The successful official response
was inspected on 2026-08-13, but a subsequent raw-response save was interrupted
by transient network refusal; no raw response hash is claimed.

The GSE63063 accession-identity audit is reproducible from
`02_ANALYSIS/audit_gse63063_overlap.py`; it uses a cached official GEO summary
with a recorded SHA-256 and file modification time. A 2026-08-13 live Entrez
recheck was blocked by the current network policy, so the status is explicitly
not upgraded to a same-day live verification.

The cohort-paper audit is reproducible from
`02_ANALYSIS/audit_source_verification.py`; its saved official PubMed responses,
record-level hashes and the title-locked GSE140829 preprint provenance are in
`03_RESULTS/source_verification_audit.*`.
