# Public MCI blood candidate-discovery extension

Access date: 2026-08-14.

## Purpose and fixed admission gate

This is an extension to `public_mci_blood_candidate_discovery_20260813.md`.
It checks whether literature-led public-data traces identify a previously
unregistered cohort that can improve the primary direct peripheral whole-blood
bulk-mRNA MCI-versus-HC transportability audit. A record is eligible only if it
has all of the following:

1. a source-resolvable public record;
2. human direct peripheral whole-blood bulk mRNA measurement;
3. usable MCI and cognitively healthy-control labels; and
4. evidence that it is not part of the AddNeuroMed development family.

Failure of any gate prevents model fitting, threshold selection, calibration,
pooling, or reporting it as a validation cohort.

## Sources queried

- NCBI PubMed E-utilities query: `"mild cognitive impairment"[Title/Abstract]
  AND "peripheral blood"[Title/Abstract] AND ("transcriptome"[Title/Abstract]
  OR "RNA-seq"[Title/Abstract] OR "gene expression profiling"[Title/Abstract])`.
  It returned 17 records on the access date.
- NCBI PubMed-to-GEO and PubMed-to-BioProject link checks for the literature
  records most likely to contain direct peripheral transcriptomics. These
  checks returned no linked GEO or BioProject records for the inspected papers.
- BioStudies/ArrayExpress free-text searches. The returned high-relevance
  human MCI records were AddNeuroMed (already registered) and GSE18309 (PBMC,
  already excluded); synonym expansion was broad and was not treated as a
  cohort admission mechanism.
- Official GEO accession viewer for `GSE165090`.
- Open full text for the two relevant literature records:
  [PMID 42529199](https://pubmed.ncbi.nlm.nih.gov/42529199/) and
  [PMID 36031604](https://pubmed.ncbi.nlm.nih.gov/36031604/).

## Literature-led candidates

### HRA019152: independent but incompatible PBMC amyloid-status context

The open report linked to [PMID 42529199](https://pubmed.ncbi.nlm.nih.gov/42529199/)
states that raw sequence data were deposited at GSA-Human accession
`HRA019152`. It describes 48 Chinese participants, six groups formed by SCD,
MCI, and dementia crossed with A-beta PET status (eight per group). The MCI
comparison is therefore eight A-beta-positive versus eight A-beta-negative
participants. Its sample preparation removes plasma and residual red cells and
isolates RNA from PBMCs, followed by poly(A)-selected bulk RNA-seq.

This is independent-source material and can potentially support a separately
preregistered MCI-internal A-beta-status exploration. It fails the primary
admission gate because it is PBMC rather than unsorted whole blood, lacks a
cognitively healthy-control group, and has only eight MCI participants per
comparison arm. Direct GSA-Human page retrieval did not yield inspectable
record metadata during this audit, so download/access status remains
`pending official-access verification`. It is not a primary or secondary
validation cohort.

### GSE165090: literature accession mismatch, not a candidate cohort

The open BMC Medicine article linked to [PMID 36031604](https://pubmed.ncbi.nlm.nih.gov/36031604/)
states that its ANM1 dataset is `GSE165090`. The official GEO accession viewer
instead identifies `GSE165090` as **Functional role of SERPINA1 in
trophoblasts**, a public trophoblast RNA-seq study. It is not an MCI blood
cohort. This independently reproducible mismatch means the article's
data-availability accession cannot establish ANM1 provenance.

This observation is not evidence of misconduct or a biological finding. It is
an example of why accession-to-paper verification must precede cohort reuse.
It adds no data to the analysis package.

## Conclusion

No new dataset meets all four fixed criteria. The newly identified HRA019152
is explicitly excluded from the whole-blood MCI-versus-HC estimand, while
GSE165090 is an unrelated accession cited in error by a secondary article.
The primary manuscript should therefore retain its reliability and
transportability framing rather than claiming that an additional positive
independent validation cohort has been found.

## Source URLs

- https://pubmed.ncbi.nlm.nih.gov/42529199/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13416258/
- https://pubmed.ncbi.nlm.nih.gov/36031604/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9422129/
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165090
- https://ngdc.cncb.ac.cn/gsa-human/
