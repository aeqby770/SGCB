# SGCB: Generalized-Gamma Differential Expression for Bulk RNA-seq

**SGCB** is an R package for two-group bulk RNA-seq differential analysis.
It fits gene-wise three-parameter Generalized Gamma (GG) models and reports
mean-shift and distribution-shift evidence in one result table.

## Key Features

- **GG model with three parameters** (`alpha`, `beta`, `gamma`) for location/scale/shape differences
- **Information-geometric fitting**: natural-gradient MAP with hierarchical prior and small-sample stabilization
- **Variance shrinkage**: TMM preprocessing + limma-style fitFDist on log2-scale pooled variance
- **Testing outputs**: primary DE (`pvalue_t`), diagnostic channels (`pvalue_mu_wald`, `pvalue_manifold`), DD channel (`pvalue_dd` from DV + DG-gamma), and omnibus score (`SGCB_Score`)
- **LFC shrinkage**: Cauchy prior (apeglm-style)
- **Contrast-defined deployment wrappers**: `sgcbContrast()` and `sgcbPairwise()` map a pre-defined two-group comparison back to the unchanged `sgcbDE()` core
- **Optional bootstrap**: calibrated confidence interval and p-value columns
- **C++/OpenMP core** via Rcpp

## Installation

```r
# Install from GitHub
devtools::install_github("aeqby770/SGCB")

# Or clone and install locally
# git clone https://github.com/aeqby770/SGCB.git
# R CMD INSTALL SGCB
```

**Requirements**: R >= 4.0.0 and Rcpp. A C++ compiler with C++17 and OpenMP support is recommended.

## Example

The example below illustrates the expected input structure and output format:

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

`sgcbDE()` is the main entry point. It performs filtering, normalization,
GG fitting, hypothesis testing, multiple-testing correction, and effect-size shrinkage.

```r
res <- sgcbDE(
  counts,              # integer count matrix (genes x samples), raw counts
  group,               # character/factor vector of length ncol(counts), exactly 2 levels
  alpha       = 0.1,   # FDR threshold used by helper print/significant filters
  use_manifold_test = TRUE,   # report manifold channel (requires n >= 6 per group)
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

## Contrast-Defined Deployment

`sgcbDE()` is the inferential core and should be used when your study is already a clean two-group comparison.

If your primary analysis starts from a richer design, but you want to follow up on one pre-defined two-group contrast, use the wrappers below:

```r
# sample_data rownames must match colnames(counts)
sample_data <- data.frame(
  condition = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
  batch = c("x", "y", "z", "x", "y", "z", "x", "y", "z"),
  row.names = colnames(counts)
)

# Follow up one target contrast already defined upstream
res_ab <- sgcbContrast(
  counts,
  sample_data = sample_data,
  group_col = "condition",
  contrast_levels = c("A", "B")
)

# Screen all pairwise contrasts within one factor
res_all <- sgcbPairwise(
  counts,
  sample_data = sample_data,
  group_col = "condition"
)
```

Use these wrappers only when:

1. a primary model or study design has already defined the biological contrast;
2. any covariate adjustment has already been handled upstream; and
3. you want SGCB as a secondary follow-up on the resulting two-group comparison.

These wrappers do **not** make SGCB a multi-factor inference engine. They only subset samples and forward the resulting two-group problem to `sgcbDE()`.

## Output Columns (46 base + 5 with bootstrap)

### Which columns should I use?

Typical usage is:

1. **Call DE genes** → filter on `padj < 0.05`
2. **Rank genes** → sort by `padj` or `SGCB_Score` (descending)
3. **Report fold changes** → use `log2FC_shrunk` (not raw `log2FoldChange`)
4. **Volcano plot** → x = `log2FC_shrunk`, y = `-log10(pvalue)`

For interpretation:

1. **Primary inferential channel** → `padj` / `pvalue_t` for DE
2. **Secondary distributional follow-up** → `dv_padj`, `dg_gamma_padj`, `padj_dd`
3. **Diagnostics only** → `pvalue_mu_wald`, `pvalue_manifold`

The DV/DG/DD channels are best treated as exploratory secondary signals. They can highlight structured non-mean changes, but they should not be treated as stand-alone biological validation without additional sensitivity analysis or replication.

### Three LFC columns — which one to use?

| Column | What it is | When to use |
|--------|-----------|-------------|
| `log2FoldChange` | Raw: mean(log2 treat) - mean(log2 ctrl) | Inspection and diagnostics |
| `log2FC_gg` | GG-model-based: log2(E[X]_treat / E[X]_ctrl) from fitted GG parameters | Use when the GG mean is the target estimand |
| **`log2FC_shrunk`** | **Cauchy-prior shrinkage of raw LFC** | **Recommended for ranking, visualization, and reporting** |

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

### Omnibus SGCB_Score (combines DE and DD channels)

| Column | Description |
|--------|-------------|
| `SGCB_Score` | -log10(Cauchy-combined p). Higher = stronger combined evidence from DE + DD |
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
| `pvalue_mu_wald` / `padj_mu_wald` | GG mean Wald channel (diagnostic; not used in SGCB_Score) |
| `pvalue_manifold` / `padj_manifold` | Manifold channel (diagnostic; NA when n < 6 per group) |
| `dv_pvalue` / `dv_padj` | Differential variability test (DV; based on log CV²) |
| `dv_log2_var_ratio` | log2(var_treat / var_ctrl) |
| `dv_log2_cv2_ratio` | log2(CV²_treat / CV²_ctrl) |
| `dv_stat` | DV test statistic |
| `dv_var_ctrl` / `dv_var_treat` | Per-group GG variance |
| `dv_cv2_ctrl` / `dv_cv2_treat` | Per-group GG CV² |
| `dg_alpha_pvalue` / `dg_alpha_padj` | Differential shape (alpha) test |
| `dg_alpha_log2_ratio` / `dg_alpha_stat` | Alpha ratio and test statistic |
| `dg_gamma_pvalue` / `dg_gamma_padj` | Differential shape (gamma) test |
| `dg_gamma_log2_ratio` / `dg_gamma_stat` | Gamma ratio and test statistic |
| `pvalue_dd` / `padj_dd` | Combined DD: Cauchy(DV + DG-gamma) |

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
# DE genes passing an adjusted p-value and effect-size threshold
de_genes <- significantGenes(res, padj_cutoff = 0.05, lfc_cutoff = 1)

# --- Distributional / variance changes ---
# Genes with significant distributional changes
dd_genes <- res[!is.na(res$padj_dd) & res$padj_dd < 0.05, ]

# Genes with significant DV signal
dv_genes <- res[!is.na(res$dv_padj) & res$dv_padj < 0.05, ]

# --- Omnibus ranking (DE + DD combined) ---
top_genes <- head(res[order(-res$SGCB_Score), ], 20)

# --- Volcano plot ---
plot(res$log2FC_shrunk, -log10(res$pvalue),
     pch = 20, cex = 0.5,
     col = ifelse(!is.na(res$padj) & res$padj < 0.05, "red", "grey60"),
     xlab = "log2 Fold Change (shrunk)", ylab = "-log10(p-value)",
     main = "SGCB Volcano Plot")
abline(v = c(-1, 1), lty = 2, col = "blue")
```

## Configuration

Most analyses can use the default configuration.

At present, the only `SGCBConfig` slot actively used by `sgcbDE()` is `bootB` (number of bootstrap replicates). Other slots (for example `maxIter`, `learningRate`, and `nCores`) are defined in the class for compatibility, but are not currently connected to the C++ core.

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

1. **Gene filtering** — remove genes with < `min_count` counts in < `min_samples` samples
2. **Normalization** — edgeR-style TMM factors + log2 workflow offset
3. **GG parameter fitting** — natural-gradient hierarchical MAP (Firth for n <= 10; gamma fixed to 1 for n <= 5)
4. **Variance shrinkage** — limma-style mean-trend fitFDist followed by bounded GG-informed prior fusion
5. **Primary DE channel** — moderated t-test (`pvalue_t` / `padj_t`)
6. **Diagnostic channels** — GG mean Wald and manifold
7. **DV/DG channels** — DV on log CV² and DG on alpha/gamma
8. **LFC shrinkage** — Cauchy prior shrinkage (`log2FC_shrunk`)
9. **Cauchy aggregation** — `pvalue_dd` from DV + DG-gamma; `SGCB_Score` from DE + DD (plus bootstrap channel when enabled)
10. *(Optional)* **Calibrated bootstrap** — additional uncertainty columns

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
│   ├── sgcb_adam_optimizer.cpp   # Legacy optimizer path
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

## Reference

If you use SGCB, please cite:

```
SGCB: Generalized-Gamma Differential Expression for Bulk RNA-seq (2025)
https://github.com/aeqby770/SGCB
```
