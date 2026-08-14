# GSE24206 Candidate Exclusion Audit

Audit date: 2026-08-12

## Decision

`GSE24206` is excluded from the primary MCI/healthy-control blood-transcriptome analysis and from external validation. It is not a compatible MCI dataset.

## Evidence checked

The locally staged GEO files were all readable to end-of-stream:

- `GSE24206_family.soft.gz`
- `GSE24206_series_matrix.txt.gz`
- `GPL570.annot.full.gz`

The official series metadata in the SOFT and series-matrix files agree on:

- Series title: `Validated Gene Expression Signatures of Idiopathic Pulmonary Fibrosis`
- Platform: `GPL570` Affymetrix Human Genome U133 Plus 2.0 Array
- Tissue: whole lung tissue
- Design: 17 IPF sample rows from 11 patients and 6 control specimens
- Repeated-sample structure: 6 patients contributed paired upper- and lower-lobe samples
- Public PubMed link recorded by GEO: `21974901`

Source URLs:

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24206
- PubMed: https://pubmed.ncbi.nlm.nih.gov/21974901/

## Compatibility gate

| Gate | Result | Reason |
|---|---|---|
| Direct peripheral blood specimen | FAIL | The series is whole-lung tissue. |
| HC/MCI/AD stage endpoint | FAIL | The labels are IPF and control, not MCI and healthy control. |
| Independent subject-level validation | FAIL | The series includes paired lung samples from the same patients. |
| Exact platform mapping | PASS, but irrelevant | GPL570 annotation is available; platform compatibility cannot rescue tissue or endpoint incompatibility. |

## Interpretation boundary

The dataset may be useful only as disease-specific or tissue-specific context. It must not be pooled into the blood MCI estimand, used for model fitting, feature selection, calibration, threshold selection, external validation, or a positive transportability claim.

This audit is a candidate-screening result. It does not change the frozen V4/V5 analyses or any existing figures and manuscript claims.
