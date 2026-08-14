# Candidate cross-stage genetic-triangulation protocol v1

Status: not executable. This is a gated design for a future public-data-only
supplementary analysis, not a preregistered result and not a workaround for the
missing independent MCI-HC transcriptome cohort.

## Fixed question

For a pre-frozen, small blood-transcriptome gene set, is there source-linked
genetic evidence that warrants reporting a cross-stage association with
Alzheimer disease? This question is distinct from MCI diagnosis, progression,
transportability and therapeutic efficacy.

## Candidate outcome

GWAS Catalog study `GCST90444373` is a public European-ancestry Alzheimer
disease study with full summary statistics available. It is a candidate outcome
only. Outcome acquisition, checksum, column audit and phenotype review must be
completed before any analysis.

## Non-negotiable admission gates

- No MCI-to-AD substitution: the outcome is labelled Alzheimer disease in every
  artifact and manuscript sentence.
- No inherited instruments or old result tables: source-linked blood cis-eQTL
  beta, SE, alleles, genome build and ancestry are required.
- No distance-only pruning: LD-reference clumping is required.
- No single-SNP causal language: a single-SNP Wald ratio may be descriptive but
  cannot carry the triangulation claim alone.
- No colocalization shortcut: each gene-locus test retains its own variants and
  inputs.
- No main-text rescue: all outcomes remain supplementary cross-stage context.

## Analysis sequence after gates pass

freeze genes -> retrieve and audit exposure -> LD clump -> harmonize ->
instrument-strength audit -> primary/sensitivity MR -> gene-locus
colocalization -> source-data and code release -> claim audit.

If any source or method gate fails, record `NO_ANALYSIS_EXECUTED` rather than
carrying forward a partial output.
