# GSE108215 Paired Intervention-Response Audit Protocol

Status: frozen before reading paired-score results on 2026-08-13.

## Purpose and scope

GSE108215 is a public randomized study of mindfulness awareness practice (MAP)
versus health education (HEP) in people with MCI. It is not an MCI-versus-HC
cohort. This protocol evaluates whether the historical frozen score has a
descriptive, within-subject molecular change over nine months. It cannot test
diagnostic transportability, clinical utility, cognitive benefit, treatment
efficacy, or causal mechanism.

## Data and units

- Source: [official GSE108215 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108215), accessed 2026-08-13.
- Assay: whole-blood Affymetrix Human Genome U133 Plus 2.0 (GPL570) processed expression matrix.
- Metadata design: 60 randomized MCI participants; the GEO narrative reports 51
  baseline and 36 post-treatment blood samples, while the official processed
  matrix contains 50 pre-treatment and 36 post-treatment columns. The audit
  uses only the 86 matrix-observed rows and treats the one-row discrepancy as
  a cohort-reporting boundary.
- Analysis unit: a subject-level post-minus-pre expression-score change. The title-derived subject token is used only inside the analysis runtime to pair rows and is not exported.
- Prespecified expected paired counts: HEP=18 and MAP=17. Any duplicate arm/subject/timepoint or parsing failure is a hard stop.

## Frozen score and preprocessing

The source panel remains exactly:

`GPI, PHF15, VAT1, EIF3B, UCP2, LAT1-3TM, C22ORF39, USP38, ZSWIM6, PRDM10, RP2, GCA`.

Official GPL570 annotation is used for exact symbols only. Case canonicalization
is permitted; aliases and substitutions are not. The expected partial mapping
is 10/12, with `PHF15` and `LAT1-3TM` absent. The full 12-gene score is
therefore blocked.

For each mapped gene, the analysis uses the median over all exact-symbol
probes of `log2(signal + 1)`. Gene-wise standardization is computed across all
86 expression rows before reading MAP/HEP or pre/post labels. MCI-peak genes
receive `+1/6` and MCI-trough genes receive `-1/6`, retaining the frozen
historical weights. No arm, timepoint, outcome, cognitive score, or follow-up
status can select probes, define weights, standardize features, or alter the
score.

## Planned summaries

The generated artifact reports only aggregate results: within-arm mean and
median paired score change, paired t-test and Wilcoxon signed-rank p-values,
an exploratory Welch comparison of MAP and HEP paired changes, and bootstrap
95% confidence intervals for each arm's mean change. These p-values are
descriptive/exploratory and no multiplicity, clinical outcome, or treatment-
effect conclusion is implied. Individual accessions, subject tokens, scores
and expression values remain outside the V5 package.

## Interpretation gate

Allowed: a label-leakage-free, partial frozen-score paired molecular-response
audit in a small public MCI intervention cohort.

Forbidden: MAP improved cognition, the intervention was effective, the score
validates MCI diagnosis, the complete 12-gene signature replicated, or any
clinical, causal, or mechanistic claim.
