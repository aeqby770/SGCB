# SGCB 0.6.0

- Reject constant positive groups using exact range checks and return unavailable DV results when no group pair has estimable positive variation. Deterministic regression tests cover these input boundaries; the statistical model and default thresholds are unchanged.
- Clarify exploratory mean-score interpretation, fitted-model posterior expected FDP versus repeated-sampling FDR, and the limits of the 0.01 sensitivity cutoff.
- Document sparse-DV eligibility, positive-count reporting, conditional interval reliability, observed-abundance CV² and suitable alternatives. Six affected unequal-depth matrices were refitted from their original streams after the eligibility correction; the R2 formulas, defaults and remaining numerical records are retained.

- At zero fitted null weight, rank arithmetic-mean effects by the continuous limit of the normalized null-mass score. The fitted effects and posterior null probabilities are retained.

- Select arithmetic mean or total relative variability with the endpoint argument of sgcbDE(). Both endpoints and gene-wise diagnostics remain the default.

