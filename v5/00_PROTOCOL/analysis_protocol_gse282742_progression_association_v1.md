# GSE282742 Exploratory Progression Association Protocol v1

## Purpose and evidence boundary

This supplementary analysis examines blood bulk RNA-seq associations with clinical progression status within public GSE282742. It is an **exploratory progression association**, not independent predictive validation and not causal evidence. It does not alter the primary AddNeuroMed MCI-versus-HC transportability result.

## Prespecified populations

The primary unit is one **subject-level representative-record analysis** row for each unique P-MCI or S-MCI participant. Public metadata has no verified visit date; therefore, the retained record is the first eligible P-MCI/S-MCI row in the source metadata order, not a biologically verified baseline visit. The expected representative-record population is 52 unique subjects: 22 P-MCI and 30 S-MCI. The analysis excludes one S-MCI row with missing age from covariate-adjusted models; this exclusion is recorded, never imputed silently. Rows labelled AD are excluded from the primary contrast.

The sensitivity population is exactly 11 verified P-MCI-to-AD pairs, matched by subject identifier. It is a **paired transition sensitivity** and not a treatment effect, causal estimate, or validation of a diagnostic/prognostic classifier.

## Primary statistical analysis

Raw expected counts are analyzed at subject level. Each sample receives a library-size counts-per-million transformation followed by `log2(CPM + 1)`. For every gene, ordinary least squares estimates the P-MCI versus S-MCI association, adjusted for continuous age and binary sex. The reported coefficient is P-MCI minus S-MCI on the log2-CPM scale, with 95% confidence interval, two-sided nominal P value, and Benjamini-Hochberg FDR across tested genes.

No gene selection, hyperparameter tuning, classifier fitting, receiver-operating characteristic analysis, or prediction performance assessment is permitted in this cohort.

## Secondary summaries

Pathway summaries are calculated from the complete adjusted representative-record gene-statistic ranking using the pinned registry in `gse282742_progression_gene_set_registry_v1.json`. They remain exploratory. Comparison with the pre-existing AddNeuroMed cross-study-stable 100-gene set is limited to signed-overlap/enrichment context. Limited or discordant overlap is reported as such and is not called replication.

For the 11 paired participants, the reported quantity is AD minus P-MCI within-subject log2-CPM change, with a paired t-test and Benjamini-Hochberg FDR. This is descriptive molecular-transition context only.

## Provenance and reporting

Every output records SHA256 hashes of the official count and subject-audit inputs, input counts, covariate handling, and the exact claim boundary. Output tables contain aggregate gene/pathway results only; subject identifiers and individual-level values are never published.
