# Manuscript narrative v2: transportability before performance

Status date: 2026-08-14

This non-destructive narrative reconstruction does not replace
`manuscript_draft_v1.md` or introduce new analysis results. Quantitative
statements and evidence boundaries are bound to `claim_evidence_ledger.md`,
`Code_Data_Figure_Mapping.md`, and `source_verification_register.md`.

## Manuscript identity

**Paper type:** Methods/reliability study with an MCI blood-transcriptome case
study.

**Reader problem:** A public-data biomarker model can retain internal
discrimination and related-study transfer while failing in a new expression
domain. Internal cross-validation cannot answer the clinical question implied by
a biomarker claim.

**Central contribution:** A leakage-aware audit that confines all fitted
decisions to development data, separates related-study transfer from independent
transport, and tests whether calibration, OOD diagnostics, relative-expression
models, conformal prediction, and cross-study stable signals restore transport.

**One-sentence argument:** In public blood-transcriptome modelling of MCI, we
show that fold-contained selection removes only modest internal optimism, whereas
the decisive reliability boundary is failure of independent transport under
severe expression-domain displacement, supported by locked sample-row and
cross-platform audits, with no evidence for a clinically validated biomarker.

**Non-negotiable boundary:** The work establishes an auditable reliability
boundary, not diagnostic utility, causal biology, therapeutic efficacy, or a
generalisable MCI signature.

## Terminology ledger

| Canonical term | First-use definition | Do not substitute with |
|---|---|---|
| MCI | mild cognitive impairment (MCI) | prodromal AD unless that endpoint is analysed |
| HC | healthy control (HC) | normal control after definition |
| AddNeuroMed development family | GSE63060/GSE63061, related study-family domains | two independent external cohorts |
| primary external transport audit | GSE140829 HC/MCI sample-row assessment | external validation; replication |
| secondary cross-platform transport audit | GSE249477 locked HC/MCI assessment | positive replication |
| sample-row audit | assessment with subject independence unresolved | subject-level validation |
| domain displacement | measured development-target expression separation | biological mechanism; batch effect |
| OOD | development-threshold out-of-distribution diagnostic | domain-membership probability |
| frozen 12-gene comparator | source-hash-locked historical comparator | newly discovered signature |
| partial score | sensitivity score missing exact symbols | complete-panel replication |
| rank/z, TSP, stable-signal audits | label-safe stress tests | domain-adaptation solution |

## Title candidates

**Recommended:** `Leakage-aware transportability audits define a reliability
boundary for public blood-transcriptome models of mild cognitive impairment`

1. `Internal discrimination does not establish transportability for public
   blood-transcriptome models of mild cognitive impairment`
2. `Domain displacement limits transport of public blood-transcriptome models
   for mild cognitive impairment`
3. `A leakage-aware audit separates study-family transfer from independent
   transport in blood-transcriptome modelling of mild cognitive impairment`

The recommended title states the object, methodological move, and bounded
consequence without implying clinical utility or an identified cause.

## Replacement abstract

Blood-transcriptome models have been proposed for molecular assessment of mild
cognitive impairment (MCI), but apparent performance can be distorted by
information leakage, related-study structure, and expression-domain shift. We
developed a public-data reliability audit that separates these failure modes. For
the primary MCI-versus-healthy-control endpoint, all fitted decisions were
confined to development data; we compared internal estimates, study-family
transfer, locked external transport, calibration, OOD diagnostics, and
label-safe feature constructions.

The nested AddNeuroMed comparator achieved an AUC of 0.795 (95% CI
0.722-0.860; n = 208), and a full-scope label-aware selection attack yielded
mean optimism of 0.0198. Bidirectional same-family transfer yielded AUCs of
0.770 and 0.730. By contrast, the locked
GSE140829 external sample-row audit gave an AUC of 0.483 (95% CI 0.420-0.542;
n = 383), Brier score 0.344, expected calibration error 0.340, and a
development-q95 OOD rate of 0.859. Development and target domains were
separable without disease labels (AUC 1.000; permutation P = 0.00020).
Rank/z, relative-expression pairs, stable-signal modelling, and conformal
prediction did not restore calibrated transport. A secondary cross-platform
audit was likewise weak and uncertain (rank/z AUC 0.576; 95% CI 0.387-0.776;
n = 41; Brier and ECE 0.488).

Leakage control alone is therefore insufficient for public biomarker claims:
internal discrimination and same-family transfer can coexist with failure under
target-domain displacement. This audit provides a reproducible framework for
locating that boundary, but does not establish a clinically validated MCI
biomarker.

## Main-results architecture

### 1. Cohort gates establish the estimand before modelling (Fig. 1; Table 1; Table S1)

**Opening claim:** To ensure that apparent sample size did not substitute for
endpoint compatibility or independence, we fixed cohort roles before model
evaluation.

Report roles only: GSE63060/GSE63061 are the related AddNeuroMed development
family; GSE140829 is the primary external *sample-row* transport audit;
GSE249477 is a secondary cross-platform stress test; GSE270454 is a small
frozen-score boundary audit; progression, methylation, intervention, and
incompatible candidates remain separate supplementary domains. Move detailed
exclusions to Table S1 and Methods.

### 2. Leakage control corrects internal optimism but does not answer transport (Fig. 2a)

**Opening claim:** To quantify one common source of internal optimism, we
contrasted fold-contained selection with a full-development label-aware
selection attack.

Report nested development AUC 0.795 (95% CI 0.722-0.860; n = 208) and mean
optimism 0.0198. Stop at the factual distinction: leakage control changed the
internal estimate modestly. Do not imply the attack measures all leakage routes.

### 3. Related-study transfer is not independent transport (Fig. 2b-c)

**Opening claim:** To distinguish study-family consistency from transportability,
we compared bidirectional AddNeuroMed transfer with the held-out GSE140829 audit.

Place same-family AUCs 0.770 and 0.730 beside GSE140829 AUC 0.483 (95% CI
0.420-0.542; n = 383). In the same paragraph, report age-only, sex-only, and
age-plus-sex external baselines (0.518, 0.452, and 0.498). Conclude narrowly:
the transcriptome comparator showed neither independent external discrimination
nor evidence of incremental discrimination over these simple baselines.

### 4. The external failure includes calibration and support mismatch (Fig. 3)

**Opening claim:** To determine whether low external discrimination co-occurred
with unreliable probabilities and development-to-target support mismatch, we
evaluated calibration, development-derived OOD distance, and selective risk.

Report Brier 0.344, ECE 0.340, and q95 OOD rate 0.859. State that OOD scores are
diagnostic rather than domain-membership probabilities. End with the conformal
stress test: mean coverage 0.669 at nominal 0.90. It is not an additional primary
endpoint.

### 5. Plausible invariant constructions fail to restore transport (Fig. 5)

**Opening claim:** To test whether platform-robust feature constructions could
recover transport without target-label fitting, we pre-specified rank/z,
relative-expression TSP, and cross-study stable-signal stress tests.

Lead with the common result: rank/z gave GSE140829 AUC 0.605 but remained
miscalibrated; TSP yielded AUC 0.548 in GSE140829 and 0.465 in GSE249477; stable
ridge yielded 0.541 and 0.524. A target-unlabeled latent CORAL stress test,
which is more permissive because it uses a batch of target expression values,
still gave GSE140829 AUC 0.544 (95% CI 0.484-0.597) and GSE249477 AUC 0.421
(95% CI 0.255-0.610). Use the GSE249477 rank/z result (AUC 0.576, 95% CI
0.387-0.776, Brier/ECE 0.488) as a separate cross-platform boundary audit.
Technical guarantees that target rows were scored before endpoint assessment
belong in Methods and Tables S12-S14 plus Table S18.

### 6. Domain and effect audits quantify the reliability boundary (Figs. 6-7)

**Opening claim:** To localize the observed transport failure without assigning
an unmeasured cause, we quantified phenotype-metadata-independent domain
separation and cross-cohort effect concordance.

First report domain-only AUC 1.000 and permutation P = 0.00020. Then report
related-family concordance (Spearman 0.673; 71.8% sign concordance) versus
GSE249477 discordance (Spearman -0.174 to -0.191; 42.6-44.4% sign concordance).
State once that these descriptive statistics cannot distinguish platform,
recruitment, sample processing, cell composition, or biology.

### 7. Score provenance and orthogonal analyses remain boundaries, not rescues (Fig. 4; Fig. S1; Tables S3, S5-S17)

**Opening claim:** To prevent historical signatures and orthogonal datasets from
being overinterpreted, we retained them as audited boundary analyses.

Report the exact 10/12 external mapping for the frozen comparator and refer to
the GSE270454 partial-score audit without calling it replication. Keep GSE208623,
progression, paired intervention context, MR, virtual perturbation, and other
computational analyses in the Supplementary Results. The main text should state
only that they provide exploratory association or hypothesis-generation context
and do not change the primary transport conclusion.

## Discussion architecture

1. **Central interpretation:** leakage-sensitive internal evaluation and
   same-family transfer can coexist with complete failure of independent
   transport; these are not conflicting results.
2. **Why the evidence is decisive:** external collapse is reinforced by poor
   calibration, high OOD burden, label-free domain separation, and failed rank/z,
   TSP, stable-signal, and conformal stress tests.
3. **What changes in practice:** cohort provenance, analysis unit, fold
   containment, assay mapping, external discrimination, calibration, and OOD
   diagnostics should be reported as one evidence ladder.
4. **What is not identified:** domain separation cannot distinguish platform,
   recruitment, sample processing, leukocyte composition, or biology; no causal
   explanation or clinical deployment follows.
5. **What resolves it:** a prespecified multi-site subject-level MCI-versus-HC
   cohort with fixed assay mapping, repeated-sample handling, blinded scoring,
   and independent clinical-utility analysis.

## Editorial QA before replacing the formal draft

- Keep `independent external transport assessment` and `sample-row audit` for
  GSE140829; never shorten either to `external validation`.
- Do not use `validation`, `replication`, `generalizable`, `causal`, or
  `clinically useful` for rank/z, TSP, stable-signal, partial-score, methylation,
  MR, or virtual-perturbation analyses.
- Keep the exact-symbol mapping boundary adjacent to each historical-panel result.
- Report GSE140829 performance intervals as sample-row intervals.
- Keep Fig. 5 visibly divided into its three stress-test modules so incompatible
  scales cannot imply a combined clinical metric.
- Preserve red `【PMID】` markers when citation matching is finalized.
- Do not claim the GEO screen exhausted all archives; it is a bounded query.

## Claim-evidence map

| Claim | Evidence | Status |
|---|---|---|
| Fold-contained selection reduces one measurable optimism route | Fig. 2a; `leakage_attack_summary.csv` | Supported, bounded to this route |
| Same-family transfer does not establish independent transport | Fig. 2b-c; `study_heldout_metrics.csv`; `metrics.csv` | Supported |
| Primary GSE140829 model did not discriminate MCI from HC | AUC 0.483, 95% CI 0.420-0.542; sample-row n = 383 | Supported, not subject-level |
| External failure co-occurred with calibration failure and OOD burden | Fig. 3; `metrics.csv`; `predictions_anonymized.csv` | Supported, diagnostic not causal |
| Alternate label-safe methods did not restore transport | Fig. 5; Tables S12-S14 | Supported, not an impossibility proof |
| Development and target domains are strongly separable | Fig. 6; AUC 1.000; P = 0.00020 | Supported, no cause assigned |
| Cross-platform effect discordance accompanies the boundary | Fig. 7; effect concordance tables | Supported, descriptive |
| The study establishes a clinical MCI biomarker | No prospective subject-level validation or utility study | Blocked |
| The study identifies a biological cause of external failure | No causal design or measured mediator | Blocked |

## Chinese editorial notes

- 主线不要写成“我们做了很多算法仍不行”，而应写成“我们用一系列不依赖目标标签的检验，确认失败发生在跨域可迁移性，而非单一泄露路径”。
- 将 GSE140829 的外部结果提到第二个结果小节，读者会更快理解后续 OOD、TSP、稳定信号和跨组学分析为何是边界诊断，而非补救性堆砌。
- GSE208623、MR、虚拟敲除和其他虚拟分析保留在补充材料，展示探索广度，但不能升格为主结论的因果证据。
- 这套叙事提高方法学、可复现性和可信 AI 期刊的匹配度；它不会把当前证据提升为 CNS 或临床验证型顶刊的证据等级。
