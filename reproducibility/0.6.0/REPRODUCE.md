# SGCB 0.6.0 reproduction

## Interpretation of the outputs

Mean `pvalue` and `padj` support exploratory ranking and screening. The posterior expected FDP bound conditions on the fitted empirical-Bayes model, its learned hyperparameters and finite catalogue. Independent mixed simulations had mean FDR 7.0%–12.2% at the primary 0.05 cutoff. The post hoc 0.01 cutoff still gave mean FDP 11.3% in the additional unequal-depth mixed count matrices; it provides no general calibration remedy. Table S28 and `results/mean_screen_thresholds.tsv` retain both operating points.

DV requires at least four nonconstant positive observations per group for eligibility. Report both group sizes and positive counts with availability, effect and interval. The fixed sparse reinforcing 10:10 scenario had 93/126 covered available intervals (73.81% at a nominal 95% level), and separate independent 10:10 matrices detected 2/11,520 true CV² changes. These results establish no universal adequate sample size. Nondetection with sparse positives or limited power leaves changes unresolved; intervals condition on the learned catalogue and local reference (Tables S29–S30).

CV² refers to observed abundance after preprocessing and normalization. Depth and sampling noise, normalization, batch and composition can contribute at fixed intrinsic biological dispersion. The package README and installed help give the same guidance, suitable uses and simpler alternatives.

This archive supplies the SGCB 0.6.0 source and Windows packages, analysis programs and independent-study results. Download `SGCB_0.6.0_study_data.zip` from the same [release](https://github.com/aeqby770/SGCB/releases/tag/v0.6.0). Extract both ZIPs into the same writable directory. The data asset creates `study_data/` with fixed inputs, saved study results and native comparator evidence; it contains no SGCB source package or installer.

## Install the current release

SGCB 0.6.0 is the installation and refitting target throughout. Start a fresh R session in the extracted directory, install the dependencies listed in the package README, then run:

```r
dir.create("library", showWarnings=FALSE)
install.packages("release/SGCB_0.6.0.zip", repos=NULL, type="win.binary", lib="library")
library(SGCB, lib.loc="library")
stopifnot(as.character(packageVersion("SGCB", lib.loc="library")) == "0.6.0")
```

The Windows binary is for R 4.5.x. On another platform, install `release/SGCB_0.6.0.tar.gz` with `type="source"` and a C++17 compiler. The independent study also requires edgeR, testthat, ggplot2, patchwork, scales and ragg; Python summaries require pandas and NumPy.

## Run the independent study

From the extracted root:

```shell
Rscript --vanilla --encoding=UTF-8 experiments/reproduce.R
```

The entry point resolves its root from its own location and loads the supplied current package from `library/`. It recreates each simulated input from its recorded replicate and random stream. Completed per-matrix and fixed-panel fit caches are skipped. The package tests and small endpoint comparisons run again. The independent runs use three threads; the six fixed panels use their recorded 16-thread setting.

- `results/final/` contains the 432 current SGCB fits, including the six eligibility-corrected depth matrices.
- `results/reference_scoring/` contains fixed comparison records for the original score convention and the same-input moment comparator. These are numerical comparison data, not another software installation. All new fits load SGCB 0.6.0.
- `experiments/null_boundary.R` evaluates both score conventions from the current package's mean model and posterior null masses. It does not load an earlier package or replace its functions at runtime.
- `results/boundary/`, `results/outer/`, `results/exposure_final/` and their TSV summaries provide boundary, background-refit and exposure-alignment results.
- `results/analysis_versions.tsv` maps each analysis to its generating record and source file. `results/r1_file_index.tsv` preserves original author build identifiers solely as record provenance.

### Fresh refit

Use a separate extraction to preserve the supplied results. Keep the fixed comparison records and omit only the completed fit directories when extracting the code archive:

```r
archive <- normalizePath("SGCB_0.6.0_reproducibility.zip", mustWork=TRUE)
fresh <- file.path(dirname(archive), "SGCB_0.6.0_clean_run")
stopifnot(!file.exists(fresh))
dir.create(fresh)
entries <- utils::unzip(archive, list=TRUE)$Name
computed <- grepl("^results/(final|protected|outer|exposure_final|boundary)/", entries)
utils::unzip(archive, files=entries[!computed], exdir=fresh)
utils::unzip("SGCB_0.6.0_study_data.zip", exdir=fresh)
setwd(fresh)
```

Install only SGCB 0.6.0 in this new directory, as above, and run `experiments/reproduce.R`. The saved comparison records provide the unchanged moment-comparator outputs and original-score reference. Fresh fits and score-convention calculations use the current implementation.

## Application analyses and figures

Figure/table regeneration reads the supplied fitted records. Model refitting is a separate operation using SGCB 0.6.0. The recorded author development build IDs select data files; they do not select an installed package.

The original saved-record default is build ID `0.5.9` for the main application/fixed-panel records and `0.5.7` for the mean-validation records. `SGCB_STUDY_RECORD_ID` and `SGCB_MEAN_RECORD_ID` can select explicitly generated replacement records. The common R setup checks the current 0.6.0 library. Saved primary fits are checked against their `content_md5` in `results/r1_file_index.tsv`; original input content and sample/group correspondence are checked independently of the author's old absolute directory.

Use one Python environment with pandas, NumPy, SciPy, matplotlib, seaborn, statsmodels, scikit-learn and openpyxl. The application R postprocessing uses data.table, matrixStats, survival, broom and ps; enrichment programs additionally use AnnotationDbi, org.Hs.eg.db, GO.db and reactome.db. SLE conditional-moment sensitivity compiles `study_data/development/57_cv_moment_regression.cpp` and requires a C++ toolchain (Rtools45 on Windows R 4.5.x). Native comparator outputs are already supplied. Separately rerunning comparator programs requires their named dependencies, including sn, MDSeq, gamlss and gamlss.dist; it is not required for saved-result figure regeneration.

From the extracted outer root, the PowerShell setup is:

```powershell
$studyRoot = (Get-Location).Path
$env:SGCB_REVISION_ROOT = Join-Path $studyRoot 'study_data'
$env:SGCB_DATA_ROOT = $env:SGCB_REVISION_ROOT
$env:SGCB_LIBRARY = Join-Path $studyRoot 'library'
$env:SGCB_THREADS = '16'
```

Use the R/Python executables from the environments in which the dependencies were installed. The following commands regenerate the saved null-partition and paired-pseudobulk summaries and the original figures/tables:

```shell
Rscript --vanilla --encoding=UTF-8 study_data/experiments/08_recount_existing_null_results.R
Rscript --vanilla --encoding=UTF-8 study_data/experiments/26_pseudobulk_paired_sensitivity.R
python -X utf8 study_data/visualization/30_publication_figures.py study_data
python -X utf8 study_data/visualization/82_new_publication_figures.py study_data
python -X utf8 study_data/experiments/29_publication_tables.py study_data
```

The three Python programs regenerate the original 14 figures and Tables S1–S27. The independent-study entry supplies Figure S9 and the TSV source data for Tables S28–S35. The formatted journal workbook is supplied separately. The KIRC postprocessing sequence is `20_kirc_patient_validation.R`, `21_external_kirc_validation.R`, `28_kirc_score_sensitivity_and_examples.R`, then `kirc_variability_score.R`; the SLE program is `sle_classical_monocytes_biology.R`. These read the supplied primary fits. Keep the system R temporary directory; the extracted directory and explicit result caches may contain spaces.

### Model refitting

New fits use the current package and are written under new result labels. For example, refit the complete SEQC input with:

```shell
Rscript --vanilla --encoding=UTF-8 study_data/experiments/45_run_current_method.R study_data/data/seqc_input.rds sgcb_0.6.0_seqc counts
```

`25_existing_nb_simulation.R` and `103_proteomics_matched_input.R` retain the original NB and protein preprocessing. After producing new fixed-panel fits, `123_evaluate_current_api.R` and `217_cv_matched_families.R` default to the current package's record ID; the latter refits Gamma/lognormal restrictions and reads its GG reference from the selected record. Set `SGCB_STUDY_RECORD_ID` explicitly when re-evaluating an existing saved fit. Preserve each dataset's recorded filtering, normalization, sample order and thread setting.

## Current-release verification

All 12 retained mean-validation inputs were independently replayed with SGCB 0.6.0 using the full original inputs and preprocessing at 16 threads: SEQC, six NB matrices, four protein spike-ins and Kang. All 72 comparisons of effects, standard errors, raw scores, BH-adjusted scores, availability and decisions were exact, covering 89,436 fitted rows. The six qPCR reference/threshold comparisons also agreed. The retained objects contribute their mean-endpoint fields to these benchmarks; their historical DV fields are outside this reuse claim.

The six fixed benchmark designs also have current-release refits matching their saved scores, DV estimates, intervals and decisions. For the retained full-model applications, direct source comparison isolates the changed mean zero-mass scoring and DV constant-positive eligibility paths. Evaluating old/new score expressions on all 31 saved posteriors gave identical scores. All 19 retained normalized inputs had no eligible-size constant positive groups. The large GTEx and disease primary fits were retained, rather than fully rerun in this verification.

An independent source build passed the 15 eligibility assertions and produced the same code/data sections and numerical kernels as the Windows package. The inherited `native_source_id` is an author development label, not a content fingerprint. Saved-record integrity uses the file-index checksums; actual source-build and numeric comparisons are supplied separately.

Independent execution from the supplied data regenerated the application figures and original table workbook, and completed the default independent-study entry with its saved fit caches. All 175 package assertions passed; 18 current endpoint comparisons had zero score, decision and missingness differences. The corresponding new local timings are verification measurements and do not replace the original paper's timings.

Evidence is in `results/retained_mean_replay_summary.tsv`, `retained_mean_all_fields_comparison.tsv`, `seqc_qpcr_replay_comparison.tsv`, `seqc_endpoint_replay_comparison.tsv`, `normal_source_path_checks.tsv`, `normal_score_actual_source_comparison.tsv`, `retained_input_eligibility.tsv`, `application_regeneration_comparison.tsv`, `independent_reproduction_validation.json` and `source_build_binary_comparison.json`. Table S34 identifies the original record IDs and the completed verification for each analysis.

## Numerical corrections and timing provenance

Constant positive groups are explicitly ineligible for DV estimation. Six original depth-control matrices (386, 387, 397, 402, 405 and 407) were refitted from their saved streams with three threads. Their DV records and derived summaries are updated; all six still have zero DV calls. For the mean endpoint, the stored scores, BH-adjusted scores and availability flags match the previous submission's records exactly. Other saved records retain their true generating metadata.

Runtime summaries preserve the original benchmark measurements on the 5700X3D. The six corrected refits used the 5900HX; their actual elapsed times remain in `seconds["SGCB"]`, and `timing_reference` preserves the measurements used for the original timing summaries. `results/eligibility_refit_comparison.tsv` records the numerical changes. Times from different hardware are kept separate.
