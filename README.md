# SGCB

SGCB estimates unconditional arithmetic-mean and relative-variability contrasts between two independent groups. Both endpoints include observed zeros.

## Installation

Start a fresh R session in the extracted 0.6.0 reproduction directory. Create a writable library and install dependencies:

```r
sgcb_library <- file.path(getwd(), "library")
dir.create(sgcb_library, showWarnings=FALSE)
.libPaths(c(sgcb_library, .libPaths()))
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "data.table",
                   "matrixStats", "mixsqp"),
                 repos="https://cloud.r-project.org", lib=sgcb_library)
```

For Windows R 4.5.x, use the supplied x86_64 binary:

```r
install.packages("release/SGCB_0.6.0.zip", repos=NULL,
                 type="win.binary", lib=sgcb_library)
```

On other platforms, or to compile locally, use the source archive and a C++17 compiler (matching Rtools on Windows):

```r
install.packages("release/SGCB_0.6.0.tar.gz", repos=NULL,
                 type="source", lib=sgcb_library)
```

Load the selected package and inspect its help:

```r
library(SGCB, lib.loc=sgcb_library)
packageVersion("SGCB")  # 0.6.0
help("sgcbDE", package="SGCB")
```

The historical public v0.5.9 assets reproduce R1. The 0.6.0 installers are supplied with `SGCB_0.6.0_reproducibility.zip`; use its local files for the current implementation.

## Two-group analysis

```r
both <- sgcbDE(counts, group, normalization="tmm")
mean_only <- sgcbDE(counts, group, normalization="tmm", endpoint="mean")
variability_only <- sgcbDE(counts, group, normalization="tmm", endpoint="variability")
```

Each single-endpoint call returns the estimates, adjusted scores and calls for that endpoint. Use the default for two endpoints, distribution-fit diagnostics and joint exploratory scores.

Supply a nonnegative gene-by-sample matrix and a two-level group factor in matching sample order. Effects compare the second factor level with the first. The groups must be independent. For prepared nonnegative abundances use `normalization="none"` and set filtering on the supplied scale. Protein missingness, biological replication and sample composition need to be handled according to the measurement design. The installed deterministic example is available with `source(system.file("examples", "minimal_workflow.R", package="SGCB"))`.

## Interpret mean-screening scores

`pvalue` is normalized cumulative posterior null mass; `padj` is its complete-family BH adjustment. The bound on posterior expected FDP conditions on the fitted empirical-Bayes model, learned hyperparameters and finite catalogue. Frequentist FDR averages realized FDP over repeated datasets. At the default 0.05 cutoff, independent mixed simulations had mean FDRs of 7.0%, 12.2% and 9.3% at 10:10, 10:50 and 50:50. Use these scores to rank and screen mean-change candidates for independent biological confirmation. `DE` and `significantGenes()` identify rows passing this exploratory screen.

The post hoc 0.01 cutoff reduced mean FDP in the primary 10:50 design, with a loss of power. In 32 additional unequal-depth mixed count matrices, mean FDP was still 11.3% and power 1.2%. Tightening the cutoff therefore provides no general frequentist calibration remedy. Table S28 and `results/mean_screen_thresholds.tsv` retain both evaluated cutoffs.

At a zero fitted global null weight, the mean score uses the continuous likelihood-ratio limit of normalized null mass. The full fit records its scoring regime. If every represented null likelihood is zero, the fit returns no ranked discoveries. This rule changes neither the exploratory interpretation nor the calibration limits above.

## Model and output reference

sgcbDE() retains zeros and estimates total-abundance CV² from positive generalized-gamma (GG) and shifted-mixture states. Strictly positive inputs use the paired positive-law posterior with nonzero frequency fixed at one. Inputs containing observed zeros fit joint positive-law/nonzero-frequency states, retain nested Gamma and lognormal support, and use sample-aligned positive-parameter and frequency influences in a matched Student reference. The mixture reference preserves the base-CV² lower bound. Paired state weights are fitted by marginal maximum likelihood and averaged with their group reversal after fitting. The current implementation conditions on those learned weights.

Regularized gene-wise GG fits provide descriptive diagnostics using a fixed N(0,1) shape prior and expected-information sampling covariance. Primary DV inference uses the paired population model.

| Output | Meaning |
|---|---|
| `log2FoldChange`, `lfcSE`, `pvalue`, `padj` | Empirical-Bayes posterior mean/SD of the log2 arithmetic-mean ratio, posterior null-mass score and full-family BH adjustment |
| `dv_log2_cv2_ratio`, `dv_se_log2_cv2`, `dv_pvalue`, `dv_padj` | Fitted-state mean log2 CV² contrast, uncertainty-scale summary, two-sided tail bound and full-family BH screening adjustment |
| `dv_ci_lower`, `dv_ci_upper` | Corresponding 95% interval for the log2 CV² contrast; observed-zero intervals invert the bounded Student reference |
| `ctrl_positive_n`, `treat_positive_n`, `ctrl_nonzero_fraction`, `treat_nonzero_fraction` | Positive counts and fractions in the default two-endpoint output |
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

The installed example uses deterministic input and writes `SGCB_minimal_results.csv`, showing effects, availability, intervals and channel membership. Population-model attributes are `population_model`, `population_posterior` and `population_timing`; gene-wise diagnostic details are in `fit_diagnostics`. The companion archive's `REPRODUCE.md` gives experiment and figure commands.

## Reproduce the analyses

The [0.6.0 reproduction archive](https://github.com/aeqby770/SGCB/releases/download/v0.6.0/SGCB_0.6.0_reproducibility.zip) contains the independent-study code and results. Six original unequal-depth matrices were refitted after the constant-positive eligibility correction, using the saved random streams. The [v0.5.9 archive](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9_reproducibility.zip) preserves the reference inputs and retained R1 analyses. Start with `REPRODUCE.md` in the 0.6.0 archive; it identifies producing versions, original timing references and the eligibility-corrected records. The repository's `reproducibility/` directory contains analysis and visualization code separately from the R package.

## Historical submitted implementation

The v0.4.1 tag preserves the originally submitted source and reproduction files. The v0.5.9 release preserves the software and example used for the first revision (R1).

The command below installs the original v0.4.1 method:

```r
remotes::install_github("aeqby770/SGCB", ref = "v0.4.1")
```

Use the SGCB 0.6.0 source or Windows package supplied here for this revision. Use v0.5.9 to reproduce R1. The v0.4.1 tag supports reproduction of the original submission.

## Input requirements and evaluation

Population DV requires at least four nonconstant positive observations per group and assumes finite second moments. This is an eligibility rule. Constant positive groups are unestimable. If no gene has estimable positive variation in both groups, the variability endpoint returns unavailable estimates for the retained genes. Reliability also depends on sample size, positive support and the effect mechanism. In the fixed sparse reinforcing 10:10 scenario, 126/200 comparisons had intervals and 93/126 covered the truth (73.81% for nominal 95% intervals). Power was 0%, 39% and 86.5% at 10:10, 10:50 and 50:50. Separate independent mixed 10:10 matrices detected 2/11,520 true CV² changes. These results establish no universal adequate sample size. Plan replication using expected detection fractions in each group.

For each gene, report both group sizes, positive counts, `dv_inference_available`, `dv_log2_cv2_ratio`, `dv_ci_lower`, `dv_ci_upper` and `dv_padj`. Single-endpoint output omits positive counts; obtain them from the corresponding input columns. Unavailable comparisons are unestimable. A nonsignificant result with sparse positives or limited power leaves variability changes unresolved. Interpret available intervals conditionally on the learned catalogue and local approximation, with the observed sparse-scenario undercoverage in mind. Tables S29–S30 provide mechanism-specific availability and coverage.

The target is CV² of observed abundance after preprocessing and normalization. Sequencing depth and sampling noise, normalization, batch effects and sample composition can contribute alongside biology. A significant DV screen alone leaves intrinsic biological heterogeneity unresolved. Under Gamma–Poisson sampling, exposure-normalized observed CV² is latent dispersion plus `1/(exposure * latent_mean)`; depth can therefore change the target at fixed biological dispersion (Table S31).

`sgcbDE()` is the current two-independent-group method. `sgcbContrast()` and `sgcbPairwise()` align metadata and call it. `sgcbDReg()`, `sgcbProtein()` and `sgcbPseudobulk()` retain the earlier distributional-regression workflow; the latter two call `sgcbDReg()`. Covariate sensitivities for the manuscript are separate conditional-moment analyses in the archived R1 code and results, rather than covariate terms in the primary two-group model.

On the fixed observed-zero inputs, v0.5.9 detected 216/1,000 true CV² changes with no finite-null or undefined-truth calls at 10:50, and 613 true changes with three false calls at 50:50. Observed-zero 10:10 yielded no discoveries. Positive-panel true/false counts were 598/1, 1,190/1 and 1,186/1 across 10:10, 10:50 and 50:50. Matched Gamma and lognormal controls report detection, effect accuracy and availability by scenario. Mean-only reference validations use the unchanged v0.5.7 mean implementation and original inputs. `REPRODUCE.md` identifies the implementation and inputs for each evaluation.

## Choose an analysis

SGCB is suited to zero-containing abundances when relative variability of the full distribution is a primary target, positive support is adequate in both groups, and across-gene distribution sharing is useful. The matched observed-zero 10:50 and 50:50 panels favor GG over Gamma and lognormal restrictions; the Gamma-generated positive alternatives favor the simpler Gamma restriction. In the independent same-input 50:50 comparison, a direct moment-based CV² test had power 76.9% and mean FDP 7.1%, versus 71.5% and 4.0% for SGCB. A direct CV² method is reasonable when moments are stable and its calibration is appropriate for the data; these results describe a power–error trade-off rather than uniform superiority.

For routine mean-focused RNA-seq, especially with covariates or systematic depth differences, use an established count workflow such as edgeR or DESeq2. Limma is a simpler mean-focused option for appropriately prepared continuous abundances. Select the method by the estimand, replication and measurement design. The two-group SGCB mean endpoint remains exploratory.
