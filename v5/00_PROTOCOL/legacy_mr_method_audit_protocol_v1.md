# Legacy cis-eQTL-to-AD MR method audit protocol v1

Status: retrospective exclusion audit. This protocol evaluates inherited output;
it does not rerun MR, generate instruments, or introduce an MCI causal result.

## Question

Can the inherited cis-eQTL-to-AD MR and colocalization output be used as causal
or MCI-stage evidence in this package?

## Required evidence gates

The inherited analysis must retain all of the following before any causal role
could be considered:

1. source-linked eQTL and outcome summary data, with the actual exposure beta
   and standard error used for each instrument;
2. LD-reference-based instrument clumping rather than distance-only pruning;
3. gene-locus-specific inputs for each colocalization calculation;
4. a phenotype compatible with the MCI-stage research question; and
5. an inspectable execution environment and log that reproduce the retained
   output.

Failure of any gate excludes the inherited output from causal interpretation.
The audit records code facts and file-provenance facts only. It does not claim
that the corresponding methods are invalid when correctly implemented.

## Online-method verification

Official TwoSampleMR documentation retrieval was attempted on 2026-08-14 but
was unavailable from the current network environment. This record is therefore
`PENDING_ONLINE_METHOD_SOURCE_VERIFICATION`; the exclusion decision is based on
the package's explicit evidence gates and the inspectable legacy code, not on
an unverified external-method claim.

## Outputs

`03_RESULTS/legacy_mr_method_audit.json` and
`03_RESULTS/legacy_mr_method_findings.csv` store code fingerprints, output
counts and the resulting claim boundary. They contain no raw genotype, eQTL or
participant-level data.
