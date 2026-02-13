# SGCB: Shallow Generalized-Gamma Calibrated Bootstrap

**SGCB** is an R package for bulk RNA-seq differential expression analysis.
It models gene expression with the three-parameter Generalized Gamma (GG)
distribution and uses information-geometric inference with empirical Bayes
variance shrinkage.

## Key Features

- **Generalized Gamma distribution**: Three-parameter family (alpha, beta, gamma) that subsumes Gamma and Weibull — captures distributional shape changes invisible to NB-based methods
- **Information-geometric inference**: Natural gradient on the Fisher manifold for efficient GG parameter estimation
- **Empirical Bayes shrinkage**: limma-style fitFDist (Smyth 2004) with mean-variance trend, plus optional GG-informed precision weighting
- **Multi-channel testing**: Moderated t-test (DE), manifold Wald chi-squared (DD), differential variance (DV), and differential shape (DG) — combined via Cauchy combination into a single SGCB_Score
- **LFC shrinkage**: Cauchy-prior shrinkage (apeglm-style)
- **Calibrated bootstrap**: Optional m-out-of-n bootstrap confidence intervals
- **OpenMP parallelization**: C++ core via Rcpp, multi-threaded

## Installation

```r
# Install from source directory
devtools::install("path/to/SGCB")

# Or from command line
# R CMD INSTALL SGCB
```

**Requirements**: R >= 4.0.0, Rcpp. A C++ compiler with C++17 and OpenMP support is recommended.

## Complete Runnable Example

Copy-paste this into R to verify your installation works:

```r
library(SGCB)

# --- 1. Simulate a count matrix (500 genes x 6 samples, 3 vs 3) ---
set.seed(42)
n_genes <- 500
counts <- matrix(rnbinom(n_genes * 6, mu = 200, size = 10),
                 nrow = n_genes, ncol = 6)
rownames(counts) <- paste0("gene_", seq_len(n_genes))
colnames(counts) <- paste0("sample_", 1:6)

# Make the first 50 genes differentially expressed (3x fold change in treat)
counts[1:50, 4:6] <- counts[1:50, 4:6] * 3

# --- 2. Define groups ---
group <- c("ctrl", "ctrl", "ctrl", "treat", "treat", "treat")

# --- 3. Run SGCB ---
res <- sgcbDE(counts, group)

# --- 4. View summary ---
print(res)

# --- 5. Extract significant genes ---
sig <- significantGenes(res, padj_cutoff = 0.05)
cat("Number of significant genes:", nrow(sig), "\n")

# --- 6. Look at top hits ---
head(sig[, c("gene_id", "log2FoldChange", "log2FC_shrunk", "pvalue", "padj",
             "SGCB_Score")])
```

## The `sgcbDE()` Function

`sgcbDE()` is the **only function you need to call**. Everything — normalization, parameter estimation, testing, multiple comparison correction, LFC shrinkage — happens inside it.

```r
res <- sgcbDE(
  counts,              # integer count matrix (genes x samples), raw counts
  group,               # character/factor vector of length ncol(counts), exactly 2 levels
  alpha       = 0.1,   # FDR threshold for significance
  use_manifold_test = TRUE,   # enable manifold distance test (recommended)
  bootstrap   = FALSE, # set TRUE for calibrated bootstrap CIs (slower)
  n_boot      = 200,   # number of bootstrap replicates (only if bootstrap=TRUE)
  min_count   = 10,    # gene filtering: minimum count in at least min_samples samples
  min_samples = 2,     # gene filtering: minimum number of samples passing min_count
  config      = NULL   # optional SGCBConfig object (NULL = default settings)
)
```

**Input**: Raw (unnormalized) count matrix. Normalization is done internally.

**Group ordering**: R's `factor()` sorts levels alphabetically. The **first level is control**, the second is treatment. For example, `c("ctrl", "treat")` → ctrl is control; but `c("WT", "KO")` → **KO becomes control** (K < W). To override, pass an explicit factor:

```r
group <- factor(c("WT","WT","WT","KO","KO","KO"), levels = c("WT", "KO"))
# Now WT = control (level 1), KO = treatment (level 2)
```

**Output**: A data.frame of class `SGCBResults` with one row per gene that passes filtering. Positive `log2FoldChange` means higher in treatment.

## Output Columns (43 base + 5 with bootstrap)

### Which columns should I use?

For most users, the workflow is:

1. **Call DE genes** → filter on `padj < 0.05`
2. **Rank genes** → sort by `padj` or `SGCB_Score` (descending)
3. **Report fold changes** → use `log2FC_shrunk` (not raw `log2FoldChange`)
4. **Volcano plot** → x = `log2FC_shrunk`, y = `-log10(pvalue)`

### Three LFC columns — which one to use?

| Column | What it is | When to use |
|--------|-----------|-------------|
| `log2FoldChange` | Raw: mean(log2 treat) - mean(log2 ctrl) | Debugging only |
| `log2FC_gg` | GG-model-based: log2(E[X]_treat / E[X]_ctrl) from fitted GG parameters | Advanced: when GG fit is trusted |
| **`log2FC_shrunk`** | **Cauchy-prior shrinkage of raw LFC** | **Use this for everything** (ranking, visualization, reporting) |

### Core DE results

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier (from rownames of input matrix) |
| `baseMean` | Mean normalized expression across both groups |
| `log2FoldChange` | Raw log2 fold change |
| `log2FC_gg` | GG-model-based log2 fold change |
| `log2FC_shrunk` | Shrunk LFC (**recommended for reporting**) |
| `lfcSE` | Standard error of the LFC |
| `stat` | Moderated t-statistic |
| `pvalue` | Primary p-value (= `pvalue_t`, moderated t-test) |
| `padj` | BH-adjusted p-value (**use this to call DE genes**) |

### Omnibus SGCB_Score (combines all test channels)

| Column | Description |
|--------|-------------|
| `SGCB_Score` | -log10(Cauchy-combined p). Higher = stronger multi-channel evidence |
| `SGCB_Score_p` | Raw omnibus p-value |
| `SGCB_Score_padj` | BH-adjusted omnibus p-value |

### Group-level summaries

| Column | Description |
|--------|-------------|
| `mean_ctrl` / `mean_treat` | Per-group normalized mean |
| `ctrl_alpha`, `ctrl_beta`, `ctrl_gamma` | Fitted GG parameters for control group |
| `treat_alpha`, `treat_beta`, `treat_gamma` | Fitted GG parameters for treatment group |
| `geodesic_dist` | Fisher-Rao geodesic distance between group GG parameters |

### Individual test channels

| Column | Description |
|--------|-------------|
| `pvalue_t` / `padj_t` | Moderated t-test (same as `pvalue`/`padj`) |
| `pvalue_mu_wald` / `padj_mu_wald` | GG-based Wald test on mean parameter |
| `pvalue_manifold` / `padj_manifold` | Manifold distance chi-squared test (NA when n < 6 per group) |
| `dv_pvalue` / `dv_padj` | Differential variance (DV) test |
| `dv_log2_var_ratio` | log2(var_treat / var_ctrl) |
| `dv_stat` | DV test statistic |
| `dv_var_ctrl` / `dv_var_treat` | Per-group GG variance |
| `dg_alpha_pvalue` / `dg_alpha_padj` | Differential shape (alpha) test |
| `dg_alpha_log2_ratio` / `dg_alpha_stat` | Alpha ratio and test statistic |
| `dg_gamma_pvalue` / `dg_gamma_padj` | Differential shape (gamma) test |
| `dg_gamma_log2_ratio` / `dg_gamma_stat` | Gamma ratio and test statistic |
| `pvalue_dd` / `padj_dd` | Combined DD: Cauchy(manifold + DV). NA when n < 6 per group |

### Bootstrap columns (only when `bootstrap = TRUE`)

| Column | Description |
|--------|-------------|
| `T_obs` | Observed LR test statistic |
| `T_null_median` | Median of null bootstrap distribution |
| `T_CI_lo` / `T_CI_hi` | 95% bootstrap CI for the test statistic |
| `pvalue_bootstrap` | Bootstrap p-value |

## How to Interpret Results

```r
# --- Standard DE analysis ---
# Significant DE genes (most common use case)
de_genes <- significantGenes(res, padj_cutoff = 0.05, lfc_cutoff = 1)

# --- Advanced: distributional / variance changes ---
# Genes with significant distributional changes (not just mean shift)
dd_genes <- res[!is.na(res$padj_dd) & res$padj_dd < 0.05, ]

# Genes with significant variance changes only
dv_genes <- res[!is.na(res$dv_padj) & res$dv_padj < 0.05, ]

# --- Omnibus ranking (all channels combined) ---
top_genes <- head(res[order(-res$SGCB_Score), ], 20)

# --- Volcano plot ---
plot(res$log2FC_shrunk, -log10(res$pvalue),
     pch = 20, cex = 0.5,
     col = ifelse(!is.na(res$padj) & res$padj < 0.05, "red", "grey60"),
     xlab = "log2 Fold Change (shrunk)", ylab = "-log10(p-value)",
     main = "SGCB Volcano Plot")
abline(v = c(-1, 1), lty = 2, col = "blue")
```

## Configuration (Advanced)

Most users can ignore this — defaults work well for typical RNA-seq.

Currently, the only SGCBConfig slot actively used by `sgcbDE()` is `bootB` (number of bootstrap replicates). Other slots (e.g., `maxIter`, `learningRate`, `nCores`) are defined in the class but are **not yet wired** to the C++ core — the C++ layer uses its own internal defaults.

```r
cfg <- SGCBConfig(bootB = 1000L)
res <- sgcbDE(counts, group, bootstrap = TRUE, config = cfg)
```

Alternatively, use the `n_boot` parameter directly (without a config object):

```r
res <- sgcbDE(counts, group, bootstrap = TRUE, n_boot = 1000)
```

See `?SGCBConfig` for the full slot list.

## What Happens Inside `sgcbDE()`

1. **Gene filtering** — removes genes with < `min_count` counts in < `min_samples` samples
2. **Normalization** — DESeq2-style median-of-ratios (Anders & Huber 2010)
3. **GG parameter fitting** — natural gradient + hierarchical Bayes MAP with Firth penalty (n <= 10) and gamma=1 reduction (n <= 5)
4. **Variance shrinkage** — limma-style fitFDist with mean-variance trend (Smyth 2004)
5. **Moderated t-test** — primary DE test (robust to model misspecification)
6. **Manifold Wald test** — chi-squared test on Fisher-Rao distance (captures distributional shifts)
7. **DV/DG tests** — differential variance and shape parameter tests
8. **LFC shrinkage** — Cauchy-prior apeglm-style shrinkage
9. **Cauchy combination** — omnibus SGCB_Score from all independent channels
10. *(Optional)* **Calibrated bootstrap** — posterior predictive parameter perturbation for CIs

## Package Structure

```
SGCB/
├── R/
│   ├── AllClasses.R       # S4 class: SGCBConfig (23 slots)
│   ├── constructors.R     # SGCBConfig(), defaultSGCBConfig()
│   ├── differential.R     # sgcbDE(), significantGenes(), print.SGCBResults()
│   ├── utils.R            # Internal: DV/DG test helpers
│   └── RcppExports.R      # Auto-generated Rcpp wrappers
├── src/                   # C++ core (OpenMP-parallelized)
│   ├── info_geometry.cpp  # Main pipeline: GG fitting + shrinkage + testing
│   ├── gg_core.cpp        # GG log-likelihood, gradient, Hessian, sampling
│   ├── bootstrap.cpp      # Resampling, BH/Holm adjustment, effect sizes
│   ├── calibrated_bootstrap.cpp  # GG MLE + LR bootstrap
│   ├── parametric_bootstrap.cpp  # Smoothed bootstrap (small samples)
│   ├── permutation.cpp    # Permutation + dropout permutation
│   ├── adaptive_inference.cpp    # Adaptive variance shrinkage
│   ├── sgcb_adam_optimizer.cpp   # Adam optimizer
│   ├── estimation_accel.cpp      # LOESS, TMM, tagwise dispersion
│   ├── matrix_ops.cpp     # Parallel matrix operations
│   ├── wide_shallow_ae.cpp       # Wide-shallow autoencoder
│   ├── conditional_ae.cpp # Conditional AE (batch correction)
│   ├── sgcb_optimized.hpp # Performance utilities
│   └── fast_special.h     # Fast digamma/trigamma/lgamma
├── tests/testthat/
├── DESCRIPTION
├── NAMESPACE
└── LICENSE
```

## Dependencies

- **Required**: R (>= 4.0.0), Rcpp, stats, parallel, methods
- **Suggested**: testthat, DESeq2, edgeR, limma

## License

MIT

## Citation

If you use SGCB, please cite:

```
SGCB: Shallow Generalized-Gamma Calibrated Bootstrap for
Bulk RNA-seq Differential Expression Analysis (2024)
