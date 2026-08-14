# Verified execution environment: V5

Status: release-candidate environment record, verified 2026-08-14.

## Scope

This record identifies the environment that executed the V5 analysis and full
QA successfully. It supports local reproducibility of the packaged public-data
audit; it does not substitute for a public repository release, immutable tag,
or archival DOI.

## Interpreter

- Runtime entry point: `E:\AI_Research_Assistant\scripts\utils\python.cmd`
- Python: `3.10.11`
- Platform used for verification: Windows

## Direct dependency pins

`requirements-v5.txt` pins every third-party package directly imported by
`02_ANALYSIS/` or `99_QA/`:

- matplotlib 3.10.8
- numpy 2.2.6
- pandas 2.3.3
- Pillow 12.2.0
- scikit-learn 1.7.2
- scipy 1.15.3
- statsmodels 0.14.6

These are direct pins rather than a claim of a portable, hash-locked binary
environment. Before an external release, export and archive the resolved
environment for the target operating system alongside the immutable source tag.

## Reproduction preflight

1. Create a clean Python 3.10 environment.
2. Install `requirements-v5.txt`.
3. Run `python -m pip check`; resolve any dependency conflicts in that clean
   environment before analysis.
4. Run `99_QA\check_runtime_environment.py` and require `status: PASS`.
5. Run the commands in the README reproduction order, ending with
   `99_QA\run_full_qa.py`.

The verified V5 run used the bundled runtime above. Its package integrity is
demonstrated by the successful full QA report and direct-dependency smoke test,
not by this document alone. A `pip check` against the bundled general-purpose
runtime currently reports unrelated optional `celloracle`/`scanpy` extras;
these packages are not imported by V5 and should not be copied into a clean
release environment.
