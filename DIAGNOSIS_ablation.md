# SGCB Ablation Study Mathematical Diagnosis Report

## Scenario: 100 genes x 3v3, NB distribution, 10% DE genes (LFC=2)

## Ablation Results
| Variant | F1 | p-value used |
|------|-----|------------|
| Full (= InfoGeom) | 0.049 | min-Bonf(p_t, p_manifold) |
| NoManifold | **0** | p_t only |
| NoDropout | 0.049 | = Full (dropout not involved in infogeom mode) |
| Fast | 0.095 | p_t(fast) + manifold |
| InfoGeom | 0.049 | = Full |

## Diagnostic Conclusions

### 1. pvalue_t completely fails under infogeom mode (F1=0)

#### Causal chain:
1. GG distribution has 3 parameters (α,β,γ); at n=3, number of parameters = sample size → **non-identifiable**
2. fit_gg_hierarchical_bayes_cpp produces extremely noisy GG parameter estimates
3. GG theoretical variance Var[X]=β²(Γ(α+2/γ)/Γ(α)-(Γ(α+1/γ)/Γ(α))²) is highly sensitive to parameters
4. adaptive_squeeze_var_gg_informed_cpp uses this unreliable theoretical variance as prior → variance shrinkage direction is wrong
5. SE too large → t-statistic too small → p-values near 1

#### Comparison with DESeq2:
- DESeq2 uses mean-variance trend (LOESS) to estimate prior variance
- Prior variance depends only on mean expression, not on per-gene distributional parameter estimates
- Therefore still has statistical power at n=3 (F1=0.024)

### 2. pvalue_manifold is the only effective component

The manifold distance test does **not** depend on variance shrinkage; it directly compares ctrl vs treat GG parameter differences.
Even if per-group GG parameter estimates are biased, the **difference** between ctrl-treat can still be detected.
Similar to how paired differences eliminate shared bias.

However, the issue is: p-value calibration of the manifold distance test relies on IQR-based χ² approximation:
- df fixed at 3
- scale estimated by Q75 / χ²(0.75, df=3)
This leads to anti-conservative p-values on null data (KS=0.365).

### 3. Dropout stability is completely absent in infogeom mode

sgcb_info_geom_inference_cpp does not internally call sde_dropout_stability_cpp.
Therefore NoDropout results are identical to Full—not because dropout is unimportant, but because it was never invoked.

### 4. Fast mode outperforms InfoGeom (F1=0.095 vs 0.049)

Fast mode (.sgcbDE_fast_impl) pvalue_t comes from sgcb_adaptive_pvalue_cpp:
- Computes variance on the **log2 scale** (not raw scale)
- Uses classic limma-style d0/s0 shrinkage (independent of GG theoretical variance)
- But also incorporates GG theoretical variance (weight_gg)

Fast mode pvalue_t provides **some** detection power; combined with manifold, F1=0.095.

## Proposed Fixes

### Priority 1: Fix pvalue_t (make it independent of GG parameter estimates)

**Option A (recommended)**: In sgcb_info_geom_inference_cpp, change variance shrinkage from
"GG-informed" to "hierarchical empirical Bayes" (similar to limma):
- Do not use GG theoretical variance as prior
- Switch to adaptive_squeeze_var_hierarchical_cpp (already implemented, uses only global log-variance distribution)
- This ensures pvalue_t works correctly even at n=3

**Option B**: When n_min <= 5, completely disable GG-informed shrinkage and fall back to classical shrinkage

### Priority 2: Fix pvalue_manifold calibration

- Do not fix df=3; use IQR to jointly estimate df and scale
- estimate_scaled_chisq_df_scale_from_iqr function exists but is not used
- Or use Genomic Control: scale = median(χ²) / χ²(0.5, df=3)

### Priority 3: Enable dropout in infogeom mode

- Add SDE dropout stability to sgcb_info_geom_inference_cpp
- Combine with t-test (Cauchy combination)

### Priority 4: log2FC estimation

- Current infogeom mode directly uses log2(mean_treat/mean_ctrl)
- No shrinkage, leading to RMSE=3.17 on SEQC data
- Should add LFC shrinkage (similar to apeglm)
