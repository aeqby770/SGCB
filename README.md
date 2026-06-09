# SGCB: Generalized-Gamma Distributional Regression for RNA-seq

**SGCB** is an R package for differential expression (DE) and differential
variability (DV) analysis of RNA-seq data. It fits gene-wise Generalized
Gamma (GG) GLMs with arbitrary design matrices and reports mean-shift and
distribution-shift evidence in one result table.

## Key Features

- **GG distributional regression** — `sgcbDReg()`: mean channel `log(β) = Xb` with arbitrary design matrix; optional dispersion channel `log(α) = Wd` for DV via likelihood-ratio test
- **Blockwise estimation**: Fisher scoring (OLS) for mean, natural gradient on Fisher–Rao manifold for shape/tail, Fisher scoring (WLS) for dispersion
- **Empirical Bayes**: limma-style fitFDist variance shrinkage + moderated t-test
- **TMM normalization** built-in (edgeR-style)
- **Covariate support**: batch effects, continuous covariates, multi-group contrasts
- **Legacy two-group API** — `sgcbDE()`: fully backward-compatible, with additional DV/DG/DD diagnostic channels
- **C++17/OpenMP core** via Rcpp — 10× faster than per-group fitting

## Installation

```r
# Install from GitHub
devtools::install_github("aeqby770/SGCB")

# Or clone and install locally
# git clone https://github.com/aeqby770/SGCB.git
# R CMD INSTALL SGCB
```

**Requirements**: R >= 4.0.0 and Rcpp. A C++ compiler with C++17 and OpenMP support is recommended.

## Quick Start — `sgcbDReg()` (recommended)

```r
library(SGCB)
set.seed(42)

# Simulate 1000 genes x 10 samples (5 vs 5)
counts <- matrix(rnbinom(1000 * 10, mu = 100, size = 5), nrow = 1000)
rownames(counts) <- paste0("gene", 1:1000)
counts[1:100, 6:10] <- counts[1:100, 6:10] * 3  # 100 DE genes (3x FC)
counts[101:150, 6:10] <- matrix(                 # 50 DV genes (higher variance)
  rnbinom(50 * 5, mu = 100, size = 1), nrow = 50)
group <- c(rep(0, 5), rep(1, 5))

# --- DE only (V1) ---
res <- sgcbDReg(counts, group = group)
head(res[res$padj < 0.05, c("gene_id", "log2FoldChange", "lfcSE", "padj")])

# --- DE + DV (V2) ---
res2 <- sgcbDReg(counts, group = group, design_disp = "auto")
head(res2[res2$padj_dv < 0.05, c("gene_id", "dv_log2ratio", "padj_dv")])

# --- With batch covariate ---
batch <- rep(c(1, 2), 5)
design <- model.matrix(~ factor(group) + factor(batch))
res3 <- sgcbDReg(counts, design = design, contrast = c(0, 1, 0))
```

## Legacy Example — `sgcbDE()` (two-group only)

```r
# sgcbDE() is the original two-group entry point (fully backward-compatible)
group <- c("ctrl", "ctrl", "ctrl", "treat", "treat", "treat")
res_legacy <- sgcbDE(counts[, 1:6], group)
sig <- significantGenes(res_legacy, padj_cutoff = 0.05)
```

## The `sgcbDReg()` Function

`sgcbDReg()` is the recommended entry point for new analyses. It supports
arbitrary design matrices, batch covariates, and optional DV testing.

```r
res <- sgcbDReg(
  counts,                # integer count matrix (genes x samples), raw counts
  design  = NULL,        # N x P design matrix (if NULL, built from group)
  group   = NULL,        # group vector (used if design is NULL)
  contrast = NULL,       # P-vector contrast for DE (default: last column)
  design_disp = NULL,    # N x Q dispersion design for DV ("auto" = same as design)
  contrast_disp = NULL,  # Q-vector contrast for DV output
  min_count   = 10,      # gene filtering threshold
  min_samples = 2        # minimum samples passing min_count
)
```

**Input**: Raw (unnormalized) count matrix. TMM normalization is done internally.

**Output**: A `data.frame` (class `SGCBDRegResults`).

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier (from rownames) |
| `baseMean` | Mean fitted expression across all samples |
| `log2FoldChange` | Contrast estimate in log2 scale |
| `lfcSE` | Standard error of log2FC |
| `stat` | Moderated t-statistic |
| `pvalue` / `padj` | DE p-value and BH-adjusted |
| `alpha` / `gamma` | Fitted GG shape and tail parameters |
| `dispersion` / `dispersion_shrunk` | Pearson dispersion (raw and EB-shrunk) |
| `pvalue_dv` / `padj_dv` | DV likelihood-ratio test (when `design_disp` given) |
| `dv_log2ratio` | log2 ratio of group-specific alpha (DV effect size) |
| `lr_stat_dv` | DV chi-square statistic |

## The `sgcbDE()` Function (legacy)

`sgcbDE()` is the original two-group entry point. It remains fully backward-compatible
and provides additional diagnostic channels (manifold distance, DG, DD, SGCB_Score).

```r
res <- sgcbDE(
  counts,              # integer count matrix (genes x samples), raw counts
  group,               # character/factor vector, exactly 2 levels
  alpha       = 0.1,   # FDR threshold for print/significant helpers
  use_manifold_test = TRUE,
  bootstrap   = FALSE,
  min_count   = 10,
  min_samples = 2
)
```

**Group ordering**: R's `factor()` sorts alphabetically. First level = control.
To override: `factor(group, levels = c("WT", "KO"))`.

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

### Output Quick Reference

| Task | Recommended columns |
|------|---------------------|
| Primary DE calling | `padj`, `pvalue_t`, `stat` |
| Effect-size reporting | `log2FC_shrunk`, `lfcSE`, `baseMean` |
| Model-based mean summary | `log2FC_gg`, `pvalue_mu_wald`, `padj_mu_wald` |
| Distributional follow-up | `dv_padj`, `dg_gamma_padj`, `padj_dd` |
| Omnibus ranking | `SGCB_Score`, `SGCB_Score_p`, `SGCB_Score_padj` |
| Parameter inspection | `ctrl_alpha`, `ctrl_beta`, `ctrl_gamma`, `treat_alpha`, `treat_beta`, `treat_gamma` |
| Diagnostic review | `pvalue_manifold`, `padj_manifold`, `geodesic_dist` |
| Bootstrap output | `T_obs`, `T_null_median`, `T_CI_lo`, `T_CI_hi`, `pvalue_bootstrap` |

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

At present, the only `SGCBConfig` slot consumed directly by `sgcbDE()` is `bootB` (number of bootstrap replicates). All other slots are stored in the object for compatibility or future extension, but are not currently read by the main two-group inference path.

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
│   ├── dreg.R             # sgcbDReg() — GG distributional regression (V1+V2)
│   ├── differential.R     # sgcbDE() — legacy two-group pipeline
│   ├── AllClasses.R       # S4 class: SGCBConfig
│   ├── constructors.R     # SGCBConfig(), defaultSGCBConfig()
│   ├── utils.R            # Internal helpers
│   └── RcppExports.R      # Auto-generated Rcpp wrappers
├── src/
│   ├── gg_dreg.cpp        # GG-DReg engine: V1 (mean GLM) + V2 (dispersion GLM + LR)
│   ├── info_geometry.cpp  # Legacy pipeline: per-group GG fit + shrinkage + testing
│   ├── gg_core.cpp        # GG log-likelihood, gradient, Hessian, Fisher info
│   ├── estimation_accel.cpp  # TMM, LOESS, tagwise dispersion
│   ├── bootstrap.cpp      # Resampling, BH/Holm adjustment, effect sizes
│   ├── fast_special.h     # Fast digamma/trigamma/lgamma
│   ├── sgcb_optimized.h   # Performance utilities
│   └── RcppExports.cpp    # Auto-generated Rcpp wrappers
├── tests/
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
Jiang H, Che A, Han Y, Han Y. SGCB: Distribution-Aware Differential Analysis
with Calibrated Mean, Variability, and Shape Channels Across Omics Modalities.
IEEE Journal of Biomedical and Health Informatics (under review), 2026.
https://github.com/aeqby770/SGCB
```
