# =============================================================================
# SGCB: S4 class definitions
# =============================================================================

# -----------------------------------------------------------------------------
# SGCBConfig: configuration class
# -----------------------------------------------------------------------------
#' @title SGCBConfig class
#' @description Stores SGCB analysis configuration parameters
#' @slot epsilon Numerical stability constant
#' @slot hiddenWidth Hidden layer width
#' @slot learningRate Learning rate
#' @slot maxIter Maximum number of iterations
#' @slot adamBeta1 Adam beta1 parameter
#' @slot adamBeta2 Adam beta2 parameter
#' @slot adamEps Adam epsilon
#' @slot weightDecay Weight decay (L2 regularization)
#' @slot gradClip Gradient clipping threshold
#' @slot lrDecayFactor Learning rate decay factor
#' @slot lrDecayPatience Patience for LR decay (epochs)
#' @slot earlyStopPatience Patience for early stopping (epochs)
#' @slot minDelta Minimum improvement for early stopping
#' @slot bootB Number of bootstrap replicates
#' @slot mGrid m search grid for calibrated bootstrap
#' @slot cbReps Calibrated bootstrap repetitions
#' @slot cbTargetCoverage Target coverage for calibrated bootstrap
#' @slot ciLevel Confidence interval level
#' @slot paramLo Lower bound for parameter clipping
#' @slot paramHi Upper bound for parameter clipping
#' @slot nDropout Number of dropout replicates
#' @slot dropoutRate Dropout rate
#' @slot nCores Number of parallel cores
#' @slot seed Random seed
#' @exportClass SGCBConfig
setClass("SGCBConfig",
  slots = c(
    epsilon = "numeric",
    hiddenWidth = "integer",
    learningRate = "numeric",
    maxIter = "integer",
    adamBeta1 = "numeric",
    adamBeta2 = "numeric",
    adamEps = "numeric",
    weightDecay = "numeric",
    gradClip = "numeric",
    lrDecayFactor = "numeric",
    lrDecayPatience = "integer",
    earlyStopPatience = "integer",
    minDelta = "numeric",
    bootB = "integer",
    mGrid = "numeric",
    cbReps = "integer",
    cbTargetCoverage = "numeric",
    ciLevel = "numeric",
    paramLo = "numeric",
    paramHi = "numeric",
    nDropout = "integer",
    dropoutRate = "numeric",
    nCores = "integer",
    seed = "integer"
  ),
  prototype = list(
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
    nCores = 1L,
    seed = 12345L
  )
)
