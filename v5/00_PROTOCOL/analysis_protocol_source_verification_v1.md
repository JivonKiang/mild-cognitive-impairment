# Cohort-Paper Source Verification Protocol

Status date: 2026-08-13.

## Purpose

Freeze a reproducible provenance record for the public cohorts used in this
reliability manuscript. Cohort-paper identity, cohort role, and primary
endpoint admissibility are checked separately.

## Evidence hierarchy

1. An official GEO accession URL identifies the public series.
2. A local GEO series-matrix header can provide an explicit `!Series_pubmed_id`.
3. A saved official PubMed E-utilities response is used to verify the linked
   PMID's title, journal, date, and DOI.
4. A record becomes `ACCESSION_AND_LINKED_PAPER_VERIFIED` only when steps 2 and
   3 agree. A title match alone is insufficient.

## Explicit non-inferences

Source linkage does not establish cohort subject independence, causal biology,
clinical utility, an externally validated biomarker, or model transportability.
`GSE140829` remains a sample-row transport audit with unresolved subject
independence until an official accession-paper link and subject map are
independently verified. Its saved official GDS-to-PubMed `elink` response has
no linked PMID as of the access date; this is a conservative non-link finding,
not evidence that no associated paper exists anywhere. `GSE249477` remains a small RNA-seq stress test and
`GSE108215` remains supplementary paired intervention context.

## Reproduction

Run `02_ANALYSIS/audit_source_verification.py`. It reads the saved official
PubMed responses in `03_RESULTS`, checks their declared titles, hashes the
available local series metadata, and writes the CSV/JSON provenance artifacts.
No network request is made during regeneration; freshness is represented by the
recorded access date rather than implied by a rerun.
