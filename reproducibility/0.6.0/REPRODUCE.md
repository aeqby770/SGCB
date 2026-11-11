# SGCB 0.6.0 reproduction

## Interpretation of the outputs

Mean `pvalue` and `padj` support exploratory ranking and screening. The posterior expected FDP bound conditions on the fitted empirical-Bayes model, its learned hyperparameters and finite catalogue. Independent mixed simulations had mean FDR 7.0%–12.2% at the primary 0.05 cutoff. The post hoc 0.01 cutoff still gave mean FDP 11.3% in the additional unequal-depth mixed count matrices; it provides no general calibration remedy. Table S28 and `results/mean_screen_thresholds.tsv` retain both operating points.

DV requires at least four nonconstant positive observations per group for eligibility. Report both group sizes and positive counts with availability, effect and interval. The fixed sparse reinforcing 10:10 scenario had 93/126 covered available intervals (73.81% at a nominal 95% level), and separate independent 10:10 matrices detected 2/11,520 true CV² changes. These results establish no universal adequate sample size. Nondetection with sparse positives or limited power leaves changes unresolved; intervals condition on the learned catalogue and local reference (Tables S29–S30).

CV² refers to observed abundance after preprocessing and normalization. Depth and sampling noise, normalization, batch and composition can contribute at fixed intrinsic biological dispersion. The package README and installed help give the same guidance, suitable uses and simpler alternatives.

This archive provides the 0.6.0 source and Windows packages, analysis code and numerical results. The public v0.5.9 release supplies the R1 reference inputs and retained analyses.

Unzip SGCB_0.6.0_reproducibility.zip into a writable directory. Tables, probe records and per-matrix results can be read directly there. To regenerate the fits, extract the contents of the public [0.5.9 reproduction archive](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9_reproducibility.zip) directly into a child directory named reference_0.5.9. Its data, results, software and release folders should then be inside that child directory. Use SGCB 0.6.0 for the independent-calibration analyses and SGCB 0.5.9 for the archived R1 analyses. Run commands from this directory. Package source is in software/SGCB, and the compiled package installs into library/SGCB. The stored design and random streams in experiments/revision_calibration.R recreate each independent simulated matrix.

## Install

Windows R 4.5.x:
```r
dir.create("library",showWarnings=FALSE)
install.packages("release/SGCB_0.6.0.zip",repos=NULL,type="win.binary",lib="library")
library(SGCB,lib.loc="library")
```

For another platform use a C++17 compiler and install release/SGCB_0.6.0.tar.gz with type="source". Package dependencies are in DESCRIPTION. To prepare the separate reference library from the extracted R1 archive, run the following from the same analysis root in R:
~~~r
dir.create("reference_0.5.9/library", recursive=TRUE, showWarnings=FALSE)
install.packages("reference_0.5.9/release/SGCB_0.5.9.tar.gz", repos=NULL, type="source", lib="reference_0.5.9/library")
~~~
On Windows R 4.5.x the [released 0.5.9 binary](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9.zip) can instead be installed to that same reference library with type="win.binary". The numerical study also uses edgeR, testthat, ggplot2, patchwork, scales and ragg, in addition to SGCB's declared dependencies. Python summaries use pandas and NumPy. The R1 numerical inputs and rerun scripts are included in the [published reproduction archive](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9_reproducibility.zip).

## Reproduce the numerical comparison

```powershell
& 'D:\R\R-4.5.3\bin\Rscript.exe' --vanilla --encoding=UTF-8 'experiments/reproduce.R'
```

The entry point resolves its root from its own location and uses the supplied libraries. It can run sequentially. The reported independent simulations used three threads per fit, and the six fixed panels used 16 threads. The replicate runner skips existing per-matrix files; use the separate-directory procedure below to regenerate fits. Table S33 reports same-input timings, which vary by machine and load.

### Fresh refit in a separate directory

For a clean rerun, start with the numerical ZIP and extract everything except the saved results into a new directory. The following R commands preserve the supplied scores for comparison.

~~~r
archive <- normalizePath("SGCB_0.6.0_reproducibility.zip", mustWork=TRUE)
fresh <- file.path(dirname(archive),"SGCB_0.6.0_clean_run")
stopifnot(!file.exists(fresh))
dir.create(fresh)
entries <- utils::unzip(archive,list=TRUE)$Name
utils::unzip(archive,files=entries[!startsWith(entries,"results/")],exdir=fresh)
setwd(fresh)
~~~

Extract the linked R1 archive into reference_0.5.9 under this fresh root and install the two package versions into their separate libraries as described above. From the fresh root, run experiments/reproduce.R. The program creates the output folders and fits the matrices using the supplied 0.6.0 implementation.

- results/final contains scores, decisions, intervals and truth fields from 432 SGCB 0.6.0 fits. Each simulated input matrix is recoverable from its replicate ID.
- results/reference_scoring contains matched scores from the R1 arithmetic-mean and DV routines. The reference package can regenerate them through the baseline argument. results/reference holds eight separate 0.5.9 package calls.
- results/boundary and results/boundary_operating_characteristics.tsv contain the unequal-depth count comparisons and 16 independent global-shift backgrounds. The global-shift transformations are paired within each background.
- results/outer is the background-refit experiment with fixed probes. It measures sensitivity to refitting cross-feature quantities. All 18 full DV refit studies use the same DV implementation as the final package.
- results/analysis_versions.tsv and results/r1_file_index.tsv map inputs, results and software versions. Comparator predictions retain their method labels. The disease-application fits are retained from R1.
- results/r1_protected_final_agreement.tsv is the per-gene agreement check on the six original fixed panels. results/final_vs_original_scoring_matrices.tsv records which of the 432 fresh matrices change and which score regime triggers the change.

Mean-screening sensitivity at cutoffs 0.05 and 0.01 is in results/mean_screen_thresholds.tsv. The summary script experiments/summarize_revision.R reads the saved adjusted scores and known truths for this calculation.

The per-probe data for Table S30 are in results/outer_refit_probes.tsv. Per-matrix exposure-alignment results for Table S31 are in results/exposure_final/.

## Figures and supplementary tables

The R numerical entry point renders Figure S9 with ggplot2 and patchwork. The six main figures and eight original supplementary figures come from the archived analyses. Tables S28–S35 report independent calibration, sparse support, refit sensitivity, known exposure, paired qPCR comparisons, endpoint costs and result provenance. `results/analysis_versions.tsv` maps the producing implementations. The supplementary workbook is a separate journal file.
