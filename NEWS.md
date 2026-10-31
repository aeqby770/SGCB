# SGCB 0.5.9

- Precompute pair-specific zero probabilities and quantile-search bounds once per CV² fit.
- Evaluate paired-population columns in batches of 256.
- Keep the statistical models, priors, filtering and screening rules unchanged.
- Update tests for unavailable DV estimates, point-null mean posteriors and control-reference mean inference.

# SGCB 0.5.8

- Fit positive-law parameters and nonzero frequency jointly for observed-zero DV inference.
- Propagate sample-aligned parameter and frequency covariance through a matched Student reference, including the mixture's base-CV² bound.
- Include nested Gamma and lognormal support during population fitting and refinement.
- Add `dv_ci_lower`, `dv_ci_upper` and score/covariance regime attributes.

# SGCB 0.5.5

- Estimate the point-at-one/Jeffreys-beta frequency prior by beta-binomial marginal maximum likelihood.
- Symmetrize fitted positive-part paired weights for CV² inference.
- Use deterministic empirical-Bayes mean inference with a sparse two-normal approximation to conditional GG location laws.
- Return posterior effects, standard deviations and null probabilities from the same fitted mean model.

# SGCB 0.5.3

- Integrate location, the arithmetic-mean correction and nonzero frequency under a joint mean model.
- Include diagonal/reverse state pairs and lognormal states in the nuisance catalogue.
- Report posterior mean effects and standard deviations, with complete-family BH adjustment of cumulative posterior null-mass scores.

# SGCB 0.5.2

- Target unconditional arithmetic-mean ratios, including zero observations.
- Allow finite-first-moment states with undefined second moments for mean inference.
- Improve mixture-weight exchanges and conditional Fourier evaluation.

# SGCB 0.5.1

- Introduce a zero mass and a positive Prentice generalized-gamma model for two independent groups.
- Estimate total-abundance CV² using paired population states and gene-specific nonzero-frequency uncertainty.
- Support prepared nonnegative abundances with `normalization = "none"`.
- Report inference availability, diagnostic fits and complete-family screening adjustments.

## Upgrading from 0.4.1

The current `sgcbDE()` method changes the primary effects and output semantics. Recompute analyses when upgrading from 0.4.1. The [README](README.md) documents installation, input requirements, output fields and the distinction between current and legacy interfaces.
