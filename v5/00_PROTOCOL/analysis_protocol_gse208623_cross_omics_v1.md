# GSE208623 Cross-Omics Concordance Protocol v1

Status: PRESPECIFIED BEFORE BETA-MATRIX ACQUISITION AND RESULT GENERATION.

## Purpose and evidence boundary

This public-data-only supplementary audit asks whether sex-adjusted MCI-minus-HC
DNA-methylation differences in the independent GSE208623 peripheral-blood-
leukocyte cohort show predeclared gene-level direction or pathway-level
concordance with the existing blood-transcriptome MCI effect audit. It is an
orthogonal, cross-omics association audit. It is not mRNA model validation,
predictive validation, diagnostic-model fitting, target-data adaptation,
pooled analysis, mediation, Mendelian randomization, cell-composition proof,
causal evidence, or therapeutic evidence.

## Data sources and acquisition gate

- GSE208623 processed beta matrix: official GEO file
  `GSE208623_Norm_beta.csv.gz`, expected byte size `363136110`.
- GSE208623 sample metadata: official GEO series-matrix header already
  source-locked in this project.
- EPIC annotation: official GPL21145 supplementary file
  `GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz`, expected byte size
  `169013044`.

The beta matrix and platform annotation are accepted only after each file has
the expected byte size, a valid gzip stream, a SHA-256 hash, and a parsable
header. The beta header must map one-to-one to all 60 official GSE208623 GSM
records. Any missing, duplicated, or unaligned sample column blocks the
analysis. Raw IDAT files are not used, because this protocol evaluates the
submitter-provided processed beta matrix and does not introduce a new
normalization pipeline.

## Population and covariates

The primary and only disease contrast is MCI minus HC: 20 MCI versus 20 HC
sample rows. AD rows are excluded before model construction. Sex is included
as a binary covariate because it is available in official metadata. Age is not
included or imputed because it is absent from the inspected official metadata.
The analysis is sample-row association only; it does not claim participant-level
clinical validation beyond the public cohort metadata.

## CpG-level association

For each CpG whose beta values are finite and lie in the closed interval [0, 1],
the analysis fits ordinary least squares on M values:

`M = log2((beta + 1e-6) / (1 - beta + 1e-6))`.

The design contains an intercept, MCI indicator, and sex indicator. The reported
MCI coefficient is MCI minus HC, with its standard error, 95% confidence
interval, two-sided P value, and Benjamini-Hochberg FDR across all estimable
CpGs. The descriptive group mean beta difference is reported separately. No
cell-composition adjustment is attempted because no source-locked reference or
prespecified deconvolution method is available for this module.

## Gene mapping and concordance

Only CpGs with one unambiguous official gene-symbol mapping are eligible for a
gene-level summary. Multi-gene, missing, or ambiguous mappings are retained in
the CpG accounting table but excluded from gene-level aggregation. A mapped
gene receives the inverse-normal weighted mean of its eligible CpG-level signed
Z statistics, where weights are one per CpG; this is a descriptive aggregation,
not an independent differential-methylation test. The gene-level methylation
direction is the sign of that aggregate.

The pre-existing transcriptome reference is the unweighted mean of the
GSE63060 and GSE63061 Hedges' g MCI-minus-HC values from the locked
cross-cohort-effect audit. Its source file `cross_cohort_effects.csv` is locked
at SHA-256 `aadeab45740f751e77f9e0294f03ddb600861f892db4415557cfabe3d7adbb41`;
the analysis fails before reading effects if the file differs. Concordance is reported only as Spearman correlation
and sign-concordance fraction in the exact mapped gene overlap, with a
two-sided label-free gene-sign permutation sensitivity. Methylation direction
is not presumed to be the inverse of expression direction, so both same-sign
and opposite-sign summaries are displayed without selecting the more favorable
orientation after results are known.

## Pathway summary

Only a source-locked, locally available pathway registry may be used. If a
registry is not source-locked, pathway enrichment remains BLOCKED; no web-tool
output or post hoc pathway database is substituted. Gene-level concordance and
CpG-level associations remain reportable independently of pathway availability.

## Outputs and reporting

Outputs contain aggregate CpG and gene statistics, source hashes, sample counts,
mapping counts, and an explicit evidence boundary. They do not publish beta
values, sample identifiers, subject tokens, model predictions, thresholds, or
individual-level data. Results are supplementary orthogonal context and cannot
upgrade the primary external transport conclusion.
