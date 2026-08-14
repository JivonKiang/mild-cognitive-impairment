# Public genetic-triangulation feasibility audit

Access date: 2026-08-14.

## Purpose

This audit asks whether current public sources support a new, reproducible
genetic-triangulation branch. It does not run MR, colocalization, TWAS, or any
genetic association test. Its role is to prevent the excluded inherited
cis-eQTL-to-AD output from being replaced by another inadequately sourced
analysis.

## Official source checks

| Source | Query or record | Result | Decision consequence |
|---|---|---|---|
| GWAS Catalog REST API 2.0 | EFO trait query `mild cognitive impairment` | Zero returned trait records | No catalogue-resolvable MCI outcome was available for a direct MCI genetic analysis. |
| GWAS Catalog REST API 2.0 | EFO `MONDO_0004975`, Alzheimer disease | Public AD studies were returned; the catalog metadata release date was 2026-08-02 | AD can only be a cross-stage outcome context, not an MCI-stage surrogate. |
| GWAS Catalog study metadata | `GCST90444373` | European-ancestry AD study, 138,558 cases and 1,013,726 controls; full summary statistics available under CC0 | Candidate AD outcome only, pending a compatible, source-linked blood cis-eQTL exposure and a preregistered analysis. |
| eQTL Catalogue API | `v3/studies` metadata endpoint | HTTP 500 on 2026-08-14 | No blood eQTL instrument source was admitted during this run. Do not substitute inherited, unverified eQTL artifacts. |

## Decision

`NO_GENETIC_ANALYSIS_EXECUTED`. The package currently has no admitted MCI GWAS
outcome and no live, source-verified blood cis-eQTL instrument source. A new MR
or colocalization result would therefore be premature.

`GCST90444373` is retained as a future candidate for a narrowly worded
cross-stage AD genetic-context analysis only. It cannot validate the MCI blood
transcriptome model, replace an independent MCI-HC cohort, identify MCI causal
biology, or rescue the primary transport result.

## Required gates before any future run

1. Freeze a small gene list before exposure or outcome association access.
2. Acquire source-linked blood cis-eQTL beta, standard error, allele and sample
   metadata from an official, accessible public source.
3. Use an LD reference matched to the exposure ancestry for clumping.
4. Harmonize alleles with palindromic-variant handling and retain every
   exclusion decision.
5. Record instrument strength, primary estimate, heterogeneity,
   leave-one-out and pleiotropy sensitivity results where the instrument count
   permits.
6. Run gene-locus-specific colocalization from retained locus inputs.
7. Label every resulting statement as cross-stage AD context unless a separate,
   compatible MCI outcome becomes available.

## Source URLs

- https://www.ebi.ac.uk/gwas/rest/api/v2/metadata
- https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits?efo_trait=mild%20cognitive%20impairment
- https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits?efo_trait=Alzheimer%20disease
- https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90444373
- http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90444001-GCST90445000/GCST90444373
- https://www.ebi.ac.uk/eqtl/api/v3/studies
