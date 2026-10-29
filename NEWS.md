# SGCB 0.5.8

Development revision, 28 October 2026.

- For observed-zero DV, fit joint positive-law/nonzero-frequency states and retain nested Gamma/lognormal seeds plus matched-log-variance counterparts.
- Propagate sample-aligned positive-parameter and frequency influences, including covariance across shared states. Match the Student reference degrees of freedom to variance-contribution heterogeneity and preserve the mixture's base-CV² bound.
- Preserve the existing strictly positive paired-population path and current control-reference arithmetic-mean model.
- Expose dv_ci_lower and dv_ci_upper and identify the score/covariance regime in result attributes.
- Refresh the six fixed benchmarks, six matched family restrictions, ten fixed null partitions and three existing disease applications. Reuse unchanged mean-only validations and competitor fits with their actual fit versions recorded.

# SGCB 0.5.5

Development revision, 22 October 2026.

- Estimate the mass of the shared point-at-one/Jeffreys-beta frequency prior by beta-binomial marginal maximum likelihood. Preserve structural all-positive observations and posterior mixture variance in both arithmetic means and CV-squared.
- Average the learned positive-part paired distribution with its group reversal before CV-squared inference. KKT diagnostics refer to the unconstrained likelihood fit before this prior symmetrization.
- Replace the default mean-channel concentration-replica MCMC with deterministic joint empirical-Bayes maximum likelihood over effect class, nuisance pair and common level. The fitted hyperparameters are treated as fixed; the previous population-uncertainty interpretation no longer applies.
- Approximate pure-GG conditional log-location laws by two normals matching four cumulants. Retain conditional first two moments for shifted-tail allocation and log-frequency uncertainty. The approximation preserves the arithmetic-mean endpoint; it does not reproduce the previous heavy-tail likelihood exactly.
- Evaluate that same normal-mixture dictionary in sparse row storage using analytic bounds on likelihood and contrast-moment omission. Keep every original gene, nuisance pair and common-level grid state in the optimization. Report full KKT and posterior consistency checks.
- Return deterministic posterior effects, standard deviations and null probabilities without obsolete Monte Carlo diagnostics. Retain exact spectral/Gibbs implementations as internal development references, outside the default mean entry point.
- Retain the finite-second-moment positive-part DV engine, its zero-aware endpoint, the inverse-CDF representation and the PocketFFT backend. The shared frequency-prior estimation changes DV posteriors and therefore requires refreshed downstream calls.
- Reuse unchanged comparator fits on the original matched inputs. Numerical storage equivalence is checked separately from statistical approximation, complete-call efficiency, empirical calibration and biological interpretation.

# SGCB 0.5.3

Updated 18 October 2026; joint arithmetic-mean inference.

- Jointly propagate b(theta), location, positivity and two class-specific population distributions, using a zero-effect class and a N(0,1) log2-effect class.
- Restore diagonal/reverse shape pairs and pure lognormal states at fitted scales, while retaining the original total prior strength. The common-level basis respects group exchange.
- DE uses cumulative posterior null mass and complete-family BH. Calibration relies on the exchangeable working population and the data-adapted finite support, without uniform fixed-nuisance validity.
- Effect and standard error are moments of that same joint posterior; log2FC_shrunk is the posterior mean without a second shrinkage step.
- Retain posterior null probabilities, batch Monte Carlo error, population diagnostics, fitted bases and RNG state. Initialize the existing project seed once; inference consumes the caller's RNG stream.
- Use Mehler/Fourier integration, exact Dirichlet prediction and lazy categorical sums. Kernel and complete-pipeline timings remain distinct. The DV model and archived comparators are retained.

# SGCB 0.5.2

Updated 17 October 2026; prepared revision source.

- DE now targets the log2 ratio of unconditional arithmetic means, including zeros. Conditional GG/shifted-GG location densities, a fitted joint nuisance distribution and positivity uncertainty replace the former moderated log-mean contrast. Existing DE effects, p-values, memberships and downstream scores require recomputation.
- Mean nuisance states require a finite first moment and permit an undefined second moment. Mean inference can use fewer than four positive observations; unavailable means remain in the full family with p=1. DE inference availability, nuisance states, timing and integration checks are explicit outputs. The offset affects descriptive log summaries.
- The DV statistical model remains the population CV-squared model described below. Its continuous Jeffreys positivity posterior differs from the point-at-one/Jeffreys mixture used by DE.
- Weight exchanges now evaluate the per-unit likelihood direction. This avoids cancellation when a tiny leaving weight makes endpoint minus current density round to zero. The original failing subproblem reaches KKT=0 in three solves and agrees exactly with the independent reference.
- Mean Fourier correlations visit only the intersection of the embedded nonzero frequency bands. On the complete 4000-gene observed-zero 10:10 input, the conditional stage measured 116.08 versus 96.99 seconds with exactly identical output fields and numerical checks. This is a stage measurement, not a complete-call speedup.
- Conditional DE tails and DV posterior tails use fitted nuisance weights; their BH screens do not establish general frequentist FDR control. Fixed inputs retain small-sample failures and scenario-specific calibration errors.

# SGCB 0.5.1

Updated 14 October 2026; local revision source, not yet publicly released.

## Changes from the submitted public 0.4.1 release

- The primary `sgcbDE()` method now analyzes two independent groups using a zero mass and a positive Prentice generalized-gamma distribution. Positive log-location, log-scale and signed shape replace the earlier bounded positive shape/scale/power fit. `Q = 0` is the lognormal limit; a fixed `N(0, 1)` shape prior is retained in the descriptive gene-wise diagnostic.
- Primary DV integrates the log2 CV-squared contrast over a learned paired population distribution of finite-moment GG and shifted-tail states, together with gene-specific beta nonzero-frequency posteriors. The reported score is a two-sided posterior-tail bound. The gene-wise penalized GG fits and covariance remain descriptive diagnostics.
- DE uses a moderated mean-log-abundance contrast independently of GG variance estimates. Its variance-prior moments now use untruncated log variances consistently with the log-F equations. Brown–Forsythe log-spread and signed-Q outputs are separate summaries; signed-Q inference is exploratory.
- Separate-channel BH uses the complete retained gene family, including unavailable tests. `DE_joint_q` and `DV_joint_q` report an exploratory adjustment across the two score types; posterior-tail DV adjustment does not imply a frequentist FDR guarantee. Default `alpha` is now 0.05 (0.1 in 0.4.1).
- `normalization = "none"` accepts prepared nonnegative abundances. Count input retains TMM normalization; the offset applies to log summaries, not the positive GG likelihood. Fit status, nonzero fractions, finite moments and inferential availability are explicit outputs. `matrixStats`, `data.table` and `mixsqp` are imports.
- Restricted mixture optimization now uses projected Newton-CG and exact Hessian-vector products, with a full paired KKT tolerance of 1e-10. Identical conditional CV-squared CDFs are reused across positive-count combinations. On the complete original 4000-feature, 100-profile observed-zero input, process CPU fell by 11.1% and peak working set by 10.4%, preserving every DE/DV call and availability flag. The original two refinement levels remain in use; historical SQP code and its MIT notice are retained.
- A deterministic installed example demonstrates the population entry point without downloads or random initialization. A complete 3:3 protein input also verifies that DE is returned when population DV is unavailable for every feature.
- Analytic moment-effect gradients avoid finite-difference steps beyond a finite-moment boundary. The optimizer reports its score at the returned iterate; mean-variance trend bins now support inputs with fewer than five features. `normalization="none"` bypasses TMM. Current API and mathematical checks cover these cases; completed gene-wise reference analyses and their fixed inputs are retained.
- The package preparation entry performs a clean local rebuild and checks loaded R/native source identity. Dependent analyses use fits from the same run. Changed stages are rerun on the affected input; unchanged fits and validated likelihood caches can be reused.

## Migration

The primary method and output semantics have changed; this is not a result-preserving upgrade. Rerun an analysis when adopting 0.5.1 and report its version. Do not combine 0.4.1 gene memberships or posterior calls with current CV-squared effects.

Use `sgcbDE(counts, group, alpha = 0.05, normalization = "tmm")` for counts, with the first factor level as reference. The old `prior_de`, `prior_dv`, `prior_dg`, `prior_sd_*`, `use_gg_variance` and `use_gamma_submodel` options are not arguments to the current entry point. Use the fields documented in README.md rather than the earlier joint posterior/ABF or DG output columns. `shape_prior_sd = 1` controls the descriptive gene-wise fit; it does not set the paired population prior.

The legacy `sgcbDReg()` and bootstrap helpers are not the revised primary DV method. Covariate adjustment and donor clustering in the accompanying manuscript use the separate joint-moment code in `development/57_cv_moment_regression.cpp` and `experiments/59_covariate_CV_sensitivity.R`, not an adjustment option in `sgcbDE()`.

## Evidence and limits

The installed primary population model yields 641 true and 15 false calls at observed-zero 50:50, with all 4,000 original features in the BH family; positive 50:50 and 10:50 yield 1,187/0 and 1,181/1. The smaller observed-zero designs and positive 10:10 remain without discoveries. These are fixed-input results, not a general calibration guarantee. The following gene-wise comparisons are retained as diagnostics of the earlier estimator.

In the matched hurdle comparison, 600 equal-CV-squared shape changes at 50 samples per group gave 6 false GG discoveries versus 36 for Gamma and 138 for lognormal. GG made no discoveries in the ten existing liver null partitions. These observations support resistance to specific shape-induced distortions, not universal FDR control.

The extra shape parameter costs accuracy and power in some settings. Overall 50:50 log2-effect RMSE was 0.910 for GG versus 0.456 for Gamma; GG made no CV-squared discoveries at 10:10. Simpler families can have greater sensitivity at retrospective, truth-dependent FDP constraints. Such curves are descriptive and cannot calibrate new data. A tested first-order prior-bias correction worsened full-scenario error and is excluded from 0.5.1.

Observed count-abundance CV-squared includes sampling noise. Sparse positive observations, heavy tails and misspecification limit moment estimation. KIRC enrichment, conditional-CV-squared analysis and reconstructed survival scores use current discovery genes. Primary survival associations remain nonsignificant in held-out TCGA and externally. Independent qPCR evidence evaluates the mean channel, not the GG contribution.

## Prepared release assets

`SGCB_0.5.1.tar.gz` contains this source, documentation, license, tests and installed example. README.md gives source installation and requirements. The archive is the prepared 0.5.1 release asset; the earlier public tag `v0.4.1` does not reproduce this revision. No public tag, release, repository submission or journal submission has been created by this preparation.
