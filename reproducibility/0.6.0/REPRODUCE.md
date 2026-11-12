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

The entry point resolves its root from its own location and loads the supplied current package from `library/`. It recreates each simulated input from its recorded replicate and random stream. Existing completed fit files are skipped. The independent runs use three threads; the six fixed panels use their recorded 16-thread setting.

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

The code archive also supplies the original application and figure programs under `study_data/experiments/`, `study_data/development/` and `study_data/visualization/`. Their software-library paths point to the current `library/` directory. For an application worker, set `SGCB_REVISION_ROOT` to the extracted `study_data` directory and `SGCB_LIBRARY` to the extracted root's `library` directory. The common setup checks SGCB 0.6.0 and the matching current source.

Original record identifiers in saved filenames are preserved so existing inputs, comparator outputs and downstream results remain traceable. They are data identifiers, rather than instructions to install a development build. Downstream application programs use the supplied saved records by default; `SGCB_STUDY_RECORD_ID` selects the record identifier when evaluating a new application refit. The original qPCR and protein mean predictions remain identified in the supplied file index.

`study_data/visualization/30_publication_figures.py`, `study_data/visualization/82_new_publication_figures.py` and `study_data/experiments/29_publication_tables.py` read the supplied study records to regenerate figures and tables. These Python programs use pandas, NumPy, SciPy, matplotlib, seaborn, statsmodels and openpyxl. `visualization/figure_S9.R` renders the independent-study figure. The journal workbook is supplied separately.

## Numerical corrections and timing provenance

Constant positive groups are explicitly ineligible for DV estimation. Six original depth-control matrices (386, 387, 397, 402, 405 and 407) were refitted from their saved streams with three threads. Their DV records and derived summaries are updated; all six still have zero DV calls. For the mean endpoint, the stored scores, BH-adjusted scores and availability flags match the previous submission's records exactly. Other saved records retain their true generating metadata.

Runtime summaries preserve the original benchmark measurements on the 5700X3D. The six corrected refits used the 5900HX; their actual elapsed times remain in `seconds["SGCB"]`, and `timing_reference` preserves the measurements used for the original timing summaries. `results/eligibility_refit_comparison.tsv` records the numerical changes. Times from different hardware are kept separate.
