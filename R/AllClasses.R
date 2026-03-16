# =============================================================================
# SGCB: S4 class definitions
# =============================================================================

# -----------------------------------------------------------------------------
# SGCBConfig: configuration class
# -----------------------------------------------------------------------------
#' @title SGCBConfig class
#' @description Stores SGCB analysis configuration parameters. In the current
#' implementation, \code{sgcbDE()} reads only the \code{bootB} slot directly.
#' The remaining slots are stored for compatibility, legacy code paths, or
#' future extensions, but are not consumed by the main two-group inference
#' routine.
#' @slot epsilon Stored numerical-stability constant
#' @slot hiddenWidth Stored hidden-layer width
#' @slot learningRate Stored learning rate
#' @slot maxIter Stored maximum number of iterations
#' @slot adamBeta1 Stored legacy optimizer parameter
#' @slot adamBeta2 Stored legacy optimizer parameter
#' @slot adamEps Stored legacy optimizer parameter
#' @slot weightDecay Stored weight-decay parameter
#' @slot gradClip Stored gradient-clipping threshold
#' @slot lrDecayFactor Stored learning-rate decay factor
#' @slot lrDecayPatience Stored patience for learning-rate decay
#' @slot earlyStopPatience Stored patience for early stopping
#' @slot minDelta Stored minimum improvement threshold
#' @slot bootB Number of bootstrap replicates used when \code{bootstrap = TRUE}
#' @slot mGrid Stored calibrated-bootstrap grid
#' @slot cbReps Stored calibrated-bootstrap repetitions
#' @slot cbTargetCoverage Stored calibrated-bootstrap target coverage
#' @slot ciLevel Stored confidence-interval level
#' @slot paramLo Stored lower bound for parameter clipping
#' @slot paramHi Stored upper bound for parameter clipping
#' @slot nDropout Stored dropout replicate count
#' @slot dropoutRate Stored dropout rate
#' @slot nCores Stored parallel-core count
#' @slot seed Stored random seed
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
