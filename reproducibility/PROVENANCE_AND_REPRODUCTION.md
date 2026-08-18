# Reproduction

Release: SGCB v0.4.1

Public files: https://github.com/aeqby770/SGCB/tree/v0.4.1
Archive: https://github.com/aeqby770/SGCB/releases/download/v0.4.1/SGCB_code_data_release_0.4.1.zip

The archive contains the source package, analysis and figure scripts, derived inputs, frozen results, `BENCHMARK_MANIFEST.tsv`, and `FINAL_DATA_TRACE.tsv`.

Entry scripts:

- SEQC/MAQC-III: `scripts/analysis/01_bulk_rnaseq/01_runner_sgcb_bulk.R`; `scripts/analysis/01_bulk_rnaseq/02_runner_other_bulk.R`.
- GTEx brain DV/null FDR: `scripts/analysis/01_bulk_rnaseq/03_runner_bulk_benchmark.R`; `scripts/analysis/04_multimodal_aggregation/03_extract_null_pvalues.R`.
- TCGA-KIRC adjusted survival: `scripts/analysis/05_tcga_kirc/01_compute_channels_and_go.R`; `scripts/analysis/05_tcga_kirc/02_dv_survival_analysis.R`; `scripts/analysis/05_tcga_kirc/03_clinical_extended_analysis.R`.

`BENCHMARK_MANIFEST.tsv` records dataset versions, filters, contrasts, null construction, seeds, and software versions. `FINAL_DATA_TRACE.tsv` maps the reported figures and tables to their frozen outputs and scripts.
