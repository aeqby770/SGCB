# =============================================================================
# SGCB: Constructor functions
# =============================================================================

# -----------------------------------------------------------------------------
# SGCBConfig constructor
# -----------------------------------------------------------------------------
#' @title Create SGCBConfig configuration object
#' @description Create SGCB analysis configuration
#' @param epsilon Numerical stability constant
#' @param hiddenWidth Hidden layer width
#' @param learningRate Learning rate
#' @param maxIter Maximum number of iterations
#' @param adamBeta1 Adam beta1 parameter
#' @param adamBeta2 Adam beta2 parameter
#' @param adamEps Adam epsilon
#' @param weightDecay Weight decay (L2 regularization)
#' @param gradClip Gradient clipping threshold
#' @param lrDecayFactor Learning rate decay factor
#' @param lrDecayPatience Patience for LR decay (epochs)
#' @param earlyStopPatience Patience for early stopping (epochs)
#' @param minDelta Minimum improvement for early stopping
#' @param bootB Number of bootstrap replicates
#' @param mGrid m search grid for calibrated bootstrap
#' @param cbReps Calibrated bootstrap repetitions
#' @param cbTargetCoverage Target coverage for calibrated bootstrap
#' @param ciLevel Confidence interval level
#' @param paramLo Lower bound for parameter clipping
#' @param paramHi Upper bound for parameter clipping
#' @param nDropout Number of dropout replicates
#' @param dropoutRate Dropout rate
#' @param nCores Number of parallel cores
#' @param seed Random seed
#' @return SGCBConfig object
#' @export
SGCBConfig <- function(
  epsilon = 1e-8,
  hiddenWidth = 5000L,
  learningRate = 0.001,
  maxIter = 500L,
  adamBeta1 = 0.9,
  adamBeta2 = 0.999,
  adamEps = 1e-8,
  weightDecay = 1e-4,
  gradClip = 5.0,
  lrDecayFactor = 0.5,
  lrDecayPatience = 20L,
  earlyStopPatience = 50L,
  minDelta = 1e-6,
  bootB = 1000L,
  mGrid = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
  cbReps = 100L,
  cbTargetCoverage = 0.5,
  ciLevel = 0.95,
  paramLo = 0.01,
  paramHi = 100.0,
  nDropout = 50L,
  dropoutRate = 0.2,
  nCores = parallel::detectCores(),
  seed = 12345L
) {

  new("SGCBConfig",
    epsilon = epsilon,
    hiddenWidth = as.integer(hiddenWidth),
    learningRate = learningRate,
    maxIter = as.integer(maxIter),
    adamBeta1 = adamBeta1,
    adamBeta2 = adamBeta2,
    adamEps = adamEps,
    weightDecay = weightDecay,
    gradClip = gradClip,
    lrDecayFactor = lrDecayFactor,
    lrDecayPatience = as.integer(lrDecayPatience),
    earlyStopPatience = as.integer(earlyStopPatience),
    minDelta = minDelta,
    bootB = as.integer(bootB),
    mGrid = mGrid,
    cbReps = as.integer(cbReps),
    cbTargetCoverage = cbTargetCoverage,
    ciLevel = ciLevel,
    paramLo = paramLo,
    paramHi = paramHi,
    nDropout = as.integer(nDropout),
    dropoutRate = dropoutRate,
    nCores = as.integer(nCores),
    seed = as.integer(seed)
  )
}

# -----------------------------------------------------------------------------
# Default configuration
# -----------------------------------------------------------------------------
#' @title Get default configuration
#' @description Return default SGCB configuration
#' @return SGCBConfig object
#' @export
defaultSGCBConfig <- function() {
  SGCBConfig()
}
