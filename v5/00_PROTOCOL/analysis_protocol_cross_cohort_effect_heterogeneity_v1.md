# Cross-Cohort Effect Heterogeneity Protocol v1

## Purpose and claim boundary

This public-data audit quantifies whether gene-level MCI-minus-HC effects are directionally stable across the available blood-transcriptome cohorts. It is a descriptive transportability and heterogeneity audit, not a diagnostic-model analysis, external predictive validation, pooled disease-effect estimate, causal analysis, or proof that technical or cell-composition differences cause observed effect differences.

## Fixed analysis universe and cohorts

Before external diagnosis labels are accessed, the eligible universe is the exact case-canonicalized gene-symbol intersection of GSE63060, GSE63061, GSE140829 and GSE249477. The four cohorts retain their own expression scales and are never expression-pooled.

The MCI-HC endpoint is analyzed independently within each cohort. GSE63060 and GSE63061 are the related AddNeuroMed development family. GSE140829 is an external sample-row audit with unresolved subject independence. GSE249477 is a small secondary independent cross-platform sample-row audit. External labels are used only after the fixed shared universe and label-free target expression matrices have been constructed; they are never used for feature selection, normalization, transformation, fitting, thresholding, calibration, or target adaptation.

## Effect and heterogeneity summaries

For each cohort and gene, the analysis reports Hedges' g for MCI minus HC, its large-sample sampling variance, and its sign. Each pairwise cohort comparison reports Spearman correlation of effects, number of shared estimable genes, and sign-concordance fraction. Development-to-external comparison uses the unweighted mean of the two AddNeuroMed effects only as a descriptive reference.

For GSE63060 and GSE63061, a two-study fixed-effect Cochran Q and I-squared summary is reported only as a descriptive within-family consistency statistic, following Higgins et al. (PMID 12958120; DOI 10.1136/bmj.327.7414.557). No pooled estimate or heterogeneity P value is reported across all four cohorts. GSE140829 is excluded from all pooled-inference calculations because subject independence is unresolved.

## Outputs and interpretation

Outputs contain only aggregate gene-level effects, pairwise summaries and source hashes. They do not export sample identifiers or expression rows. The permitted conclusion is limited to the observed degree of cross-cohort effect concordance or discordance within the mapped shared universe. The audit cannot establish biological mechanism, assay causality, blood-cell causality, or clinical utility.
