# GSE127711 Candidate Exclusion Audit

Audit date: 2026-08-12

## Decision

`GSE127711` is excluded from the primary MCI-versus-HC blood-transcriptome endpoint. It may be retained only as orthogonal context for memory-related blood biology or risk phenotypes.

## Evidence checked

The NCBI GEO official SOFT metadata were read from the GSE127711 family record. GEO identifies:

- Series title: `Blood Biomarkers for Memory: Towards Early Detection of Risk for Alzheimer Disease, Pharmacogenomics, and Repurposed Drugs`
- Platform: `GPL570`
- Sample taxon: human
- Public PubMed record: `31792364`
- Study design: longitudinal within-subject blood expression studies tracking short-term memory retention measured by the Hopkins Verbal Learning Test, followed by testing in independent cohorts
- Diagnostic context: the SOFT design lists schizophrenia, schizoaffective disorder, bipolar disorder, mood, psychotic disorder NOS, major depressive disorder and PTSD

The PubMed record was independently checked through NCBI E-utilities:

- https://pubmed.ncbi.nlm.nih.gov/31792364/
- DOI: `10.1038/s41380-019-0602-2`

## Compatibility gate

| Gate | Result | Reason |
|---|---|---|
| Blood mRNA assay | PASS | GPL570 blood expression data are present. |
| Direct HC/MCI stage endpoint | FAIL | The primary phenotype is quantitative HVLT memory retention and risk-related testing, not a clean MCI/HC label. |
| Independent analysis unit | FAIL for primary use | The study explicitly includes longitudinal within-subject designs and independent cohorts with different roles. |
| Confounding by psychiatric diagnosis | FAIL for primary use | Multiple psychiatric diagnostic groups are part of the study design. |

## Interpretation boundary

GSE127711 must not be pooled with AddNeuroMed or GSE140829, used for MCI/HC model fitting, thresholding, calibration, or positive external validation. A future supplementary analysis could ask whether an MCI-derived score associates with the quantitative HVLT phenotype, but that would be a different estimand and would require a prespecified subject-level design with explicit covariate handling.
