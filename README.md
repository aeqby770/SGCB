# SGCB

Version 0.5.9 analyzes unconditional arithmetic means and relative variability in two independent groups.

The current method is available through `sgcbDE()`. [Release v0.5.9](https://github.com/aeqby770/SGCB/releases/tag/v0.5.9) provides the source archive and the Windows R 4.5.x binary. The originally submitted implementation is retained under the historical GitHub tag v0.4.1.

## Install version 0.5.9

Start a fresh R session. Set up a writable project library and install the dependencies. Windows uses CRAN binaries for these dependencies:

```r
sgcb_library <- file.path(getwd(), "library")
dir.create(sgcb_library, showWarnings = FALSE)
.libPaths(c(sgcb_library, .libPaths()))

install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "matrixStats",
                   "data.table", "mixsqp"),
                 repos = "https://cloud.r-project.org", lib = sgcb_library,
                 type = ifelse(.Platform$OS.type == "windows", "binary",
                               getOption("pkgType")))
```

Choose one of the following archive installations.

**Windows R 4.5.x, precompiled package.** [SGCB_0.5.9.zip](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9.zip) was built with R 4.5.3 for x86_64 Windows. Install directly from the release URL; compilation is unnecessary:

```r
install.packages(
  "https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9.zip",
  repos = NULL, type = "win.binary", lib = sgcb_library
)
```

**Source installation.** Install the fixed GitHub tag on other platforms or to compile locally. This requires a C++17 compiler; Windows users need Rtools matching their R version. The verified Windows toolchain is R 4.5.3 with Rtools45. See [R's source-installation instructions](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/install.packages.html).

```r
install.packages("remotes", repos = "https://cloud.r-project.org",
                 lib = sgcb_library)
remotes::install_github("aeqby770/SGCB", ref = "v0.5.9",
                        lib = sgcb_library, dependencies = NA,
                        upgrade = "never", build_vignettes = FALSE)
```

The release also supplies [SGCB_0.5.9.tar.gz](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9.tar.gz). After downloading it and installing the dependencies above, use `install.packages("SGCB_0.5.9.tar.gz", repos = NULL, type = "source", lib = sgcb_library)`.

Load the installed package and confirm its version:

```r
library(SGCB, lib.loc = sgcb_library)
packageVersion("SGCB")  # 0.5.9
find.package("SGCB")
```

## Run the installed example

```r
source(system.file("examples", "minimal_workflow.R", package = "SGCB"))
```

The example analyzes a fixed 120-gene, 100-sample matrix and writes `SGCB_minimal_results.csv`. It demonstrates the public API and its output fields. The input is deterministic and requires no download.

For your own data, provide a nonnegative gene-by-sample matrix and a two-level group factor in matching sample order. The first level is the control; effects compare the second level with the first.

```r
result <- sgcbDE(counts, group, normalization = "tmm")
result[, c("gene_id", "log2FoldChange", "padj", "dv_log2_cv2_ratio",
           "dv_padj", "dv_ci_lower", "dv_ci_upper", "dv_inference_available",
           "DE", "DV", "channel")]
de_genes <- result[result$DE, ]
dv_genes <- result[result$DV, ]
```

For prepared nonnegative abundances, use `normalization = "none"` and set filtering on their supplied scale. Positive intensities may use `offset = 0`; retain a positive offset for descriptive log summaries when the input includes zeros.

## Choose the entry point

| Entry point | Role |
|---|---|
| `sgcbDE()` | Current two-independent-group arithmetic-mean and total-abundance CV² analysis |
| `sgcbContrast()`, `sgcbPairwise()` | Metadata-aligned wrappers around the current `sgcbDE()` method |
| `sgcbDReg()`, `sgcbProtein()`, `sgcbPseudobulk()` | Retained distributional-regression implementation; the protein and pseudobulk helpers call `sgcbDReg()` |
| `SGCBConfig()` and bootstrap helpers | Configuration and utilities retained for the earlier implementation |

Use `sgcbDE()` for the v0.5.9 analysis, including independently replicated pseudobulk counts and prepared protein abundances. Aggregation, normalization, protein missingness and covariate handling follow the data's measurement and replication design. The two-group entry point requires independent groups. The companion reproduction archive supplies the manuscript's separate conditional-CV² sensitivity for covariates and donor clustering.

Version 0.5.9 changes the inference and output meanings relative to v0.4.1. Recompute an analysis when adopting the current method; [NEWS.md](NEWS.md) describes the argument and output changes.

## Model and output reference

sgcbDE() retains zeros and estimates total-abundance CV² from positive generalized-gamma (GG) and shifted-mixture states. Strictly positive inputs use the paired positive-law posterior with nonzero frequency fixed at one. Inputs containing observed zeros fit joint positive-law/nonzero-frequency states, retain nested Gamma and lognormal support, and use sample-aligned positive-parameter and frequency influences in a matched Student reference. The mixture reference preserves the base-CV² lower bound. Paired state weights are fitted by marginal maximum likelihood and averaged with their group reversal after fitting. The current implementation conditions on those learned weights.

Regularized gene-wise GG fits provide descriptive diagnostics using a fixed N(0,1) shape prior and expected-information sampling covariance. Primary DV inference uses the paired population model.

| Output | Meaning |
|---|---|
| `log2FoldChange`, `lfcSE`, `pvalue`, `padj` | Empirical-Bayes posterior mean/SD of the log2 arithmetic-mean ratio, posterior null-mass score and full-family BH adjustment |
| `dv_log2_cv2_ratio`, `dv_se_log2_cv2`, `dv_pvalue`, `dv_padj` | Fitted-state mean log2 CV² contrast, uncertainty-scale summary, two-sided tail bound and full-family BH screening adjustment |
| `dv_ci_lower`, `dv_ci_upper` | Corresponding 95% interval for the log2 CV² contrast; observed-zero intervals invert the bounded Student reference |
| `ctrl_nonzero_fraction`, `treat_nonzero_fraction` | Observed positive fractions |
| `ctrl_mu`, `ctrl_sigma`, `ctrl_Q` and treatment counterparts | Descriptive gene-wise positive-part Prentice GG parameters |
| `fit_converged`, `fit_covariance_available`, `dv_moments_finite`, `dv_inference_available` | Separate gene-wise fit diagnostics and population CV² availability |
| `logspread_pvalue`, `logspread_padj` | Separate median-centered Brown–Forsythe log-spread summary |
| `shape_Q_difference`, `shape_Q_pvalue` | Exploratory signed-Q contrast |
| `DE_joint_q`, `DV_joint_q` | Exploratory BH adjustment across the posterior-null-mass DE and posterior-tail DV scores |

The DE model uses pure finite-first-moment GG states and integrates the state-specific arithmetic-mean correction, location and point-at-one/Jeffreys-beta frequency uncertainty. A shared prior describes the control-group arithmetic mean and ordered GG state pair under both hypotheses. The global point-null fraction and learned Gaussian nonzero-effect law are independent of that baseline. Pure-GG conditional location laws use two normal components matching four cumulants; log-frequency uncertainty enters through its first two moments. A positive observation in each group is required for a DE score. Complete-family BH is applied to cumulative fitted posterior null mass. Posterior summaries condition on the learned catalogue, weights and effect law. The default branches use deterministic calculations.

The mean dictionary uses sparse row storage with an analytic likelihood/moment omission budget. `mean_integration_checks` reports the KKT gap, storage size, posterior identity check and posterior-null omission bound. This storage approximation is separate from the conditional-location model approximation.

Group levels follow `factor(group)`; effects compare the second level with the first. Specify an ordered factor to control direction. Separate-channel BH defines `DE`, `DV`, and `channel`. Interpret category labels alongside `dv_inference_available`, which identifies eligible DV estimates. With observed zeros, DV scores use the matched local reference; strictly positive scores use the positive-law posterior. Their calibration is assessed on the supplied fixed controls and null partitions. Combined-family adjustment is exploratory.

Count defaults are `normalization="tmm"`, `offset=0.5`, `min_count=10`, `min_samples=2`, `alpha=0.05`, and `shape_prior_sd=1`. The offset is used only for descriptive log-scale and log-spread summaries. Prepared nonnegative abundances use `normalization="none"`; positive intensities can use `offset=0`, while data containing zeros require a positive offset for the log summaries. Handle protein missingness and biological replication before fitting.

DE state pairs, posterior summaries, numerical checks and component times are recorded in `mean_model`, `mean_posterior`, `mean_integration_checks` and `mean_timing`.

The installed example uses deterministic input and writes `SGCB_minimal_results.csv`, showing effects, availability, intervals and channel membership. Population-model attributes are `population_model`, `population_posterior` and `population_timing`; gene-wise diagnostic details are in `fit_diagnostics`. The companion archive's `REPRODUCE.md` gives experiment, figure and document commands.

## Historical submitted implementation

The repository's `reproducibility/` directory records the original submission. The v0.4.1 tag preserves that implementation; the current method and installed example are versioned under v0.5.9.

The command below installs the original v0.4.1 method:

```r
remotes::install_github("aeqby770/SGCB", ref = "v0.4.1")
```

Use the v0.5.9 source or Windows archive for the current method. The v0.4.1 tag supports reproduction of the original submission.

## Input requirements and evaluation

Population DV requires at least four positive observations per group and assumes finite second moments. When every gene lacks eligible DV observations, the function returns DE estimates with DV marked unavailable. Sparse positives, influential tails and model misspecification affect power and calibration. The target is variability in observed abundance. `sgcbDE()` fits two independent groups; the companion reproduction archive supplies `development/57_cv_moment_regression.cpp` and the executable example `experiments/sle_classical_monocytes_biology.R` for conditional-moment CV² sensitivity with covariates and donor clustering. The KIRC and liver sensitivity fits and sample designs are also included in that archive. The interface table above identifies the legacy regression and bootstrap helpers.

On the fixed observed-zero inputs, v0.5.9 detected 216/1,000 true CV² changes with no finite-null or undefined-truth calls at 10:50, and 613 true changes with three false calls at 50:50. Observed-zero 10:10 yielded no discoveries. Positive-panel true/false counts were 598/1, 1,190/1 and 1,186/1 across 10:10, 10:50 and 50:50. Matched Gamma and lognormal controls report detection, effect accuracy and availability by scenario. Mean-only reference validations use the unchanged v0.5.7 mean implementation and original inputs. `REPRODUCE.md` identifies the implementation and inputs for each evaluation.
