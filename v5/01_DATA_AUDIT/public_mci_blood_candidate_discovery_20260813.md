# Public GEO MCI blood candidate discovery

Access date: 2026-08-13.

## Query scope and saved evidence

This is a reproducible screen of one explicit public NCBI GEO query, not a
claim that every public archive has been exhausted. The exact query was:

```text
"mild cognitive impairment" AND blood AND GSE
```

The saved query returned 38 unique GEO series metadata records. The raw NCBI
E-utilities artifacts are retained in:

- `01_DATA_AUDIT/gds_mild_cognitive_impairment_blood_search_20260813.json`
- `01_DATA_AUDIT/gds_mild_cognitive_impairment_blood_summary_20260813.json`

Narrower whole- or peripheral-blood metadata cross-checks are retained in:

- `01_DATA_AUDIT/gds_mci_whole_or_peripheral_blood_search_20260813.json`
- `01_DATA_AUDIT/gds_mci_whole_or_peripheral_blood_summary_20260813.json`

Eighteen of the 38 records were absent from the pre-existing candidate registry.
Each was screened for human specimen, direct peripheral whole-blood bulk mRNA
assay, usable MCI and HC comparator, and cohort independence from the
AddNeuroMed development family. No previously unregistered record met all of
these prespecified criteria.

## Targeted terminology expansion

A second official GEO screen was performed to reduce dependence on the literal
phrase `"mild cognitive impairment"`. The exact query was:

```text
((amnestic MCI OR aMCI OR "cognitive impairment") AND (whole blood OR peripheral blood) AND Homo sapiens) AND GSE
```

The search returned 25 GSE metadata records. Its saved official E-utilities
artifacts are:

- `01_DATA_AUDIT/gds_targeted_amci_peripheral_blood_search_20260813.json`
- `01_DATA_AUDIT/gds_targeted_amci_peripheral_blood_summary_20260813.json`

The only newly identified direct blood cohort with explicit MCI and cognitively
healthy-control labels was GSE208623 (20 MCI, 20 control, 20 AD). It is an
Illumina EPIC DNA-methylation array, not bulk mRNA expression, and is therefore
excluded from the primary mRNA estimand and all model fitting, tuning, pooling,
or performance reporting. The remaining records were already registered,
non-mRNA, CSF/PBMC/purified-cell, non-MCI, or otherwise incompatible. This
targeted screen likewise found no additional independent whole-blood bulk-mRNA
MCI-versus-HC validation cohort.

## Compatible records and fixed roles

| Record(s) | Role in this package | Why this role is limited |
|---|---|---|
| GSE63060/GSE63061 | AddNeuroMed development family | Related study family; development only |
| GSE140829 | Sample-row transport audit | Direct whole-blood microarray MCI/HC comparison, but subject independence is unresolved |
| GSE249477 | Secondary independent cross-platform transport audit | Direct blood RNA-seq MCI/HC comparison; small sample-row assessment and platform shift; not external validation or replication |
| GSE270454 | Small independent frozen-score boundary audit | Direct whole-blood RNA-seq, but only 10/12 exact frozen symbols and small MCI subgroup; no primary-model admission |

The nearby but ineligible records were retained in
`00_PROTOCOL/candidate_compatibility_registry.csv` with their assay, specimen,
endpoint, study-family, overlap, or label limitations. PBMC, purified-cell,
CSF, plasma/EV, miRNA, single-cell, intervention, progression and repeated-row
studies were not pooled into the primary direct peripheral whole-blood bulk mRNA
MCI-versus-HC estimand.

## Implication

This saved-GEO-query screen did not identify an additional independent
whole-blood bulk-mRNA MCI-versus-HC validation cohort. It strengthens the
cohort-accounting boundary, but it does not convert GSE249477 or any
incompatible record into a positive replication cohort.
