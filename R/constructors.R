# =============================================================================
# SGCB: Constructor functions
# =============================================================================

# -----------------------------------------------------------------------------
# SGCBConfig constructor
# -----------------------------------------------------------------------------
#' @title Create SGCBConfig configuration object
#' @description Create an \code{SGCBConfig} object. In the current \code{sgcbDE}
#' pipeline, only \code{bootB} is consumed directly by the analysis. All other
#' slots are retained for object compatibility, future extension, or legacy code
#' paths, but are not currently read by the main two-group inference routine.
#' @param epsilon Stored numerical-stability constant (currently not used by \code{sgcbDE})
#' @param hiddenWidth Stored hidden-layer width (currently not used by \code{sgcbDE})
#' @param learningRate Stored learning rate (currently not used by \code{sgcbDE})
#' @param maxIter Stored iteration limit (currently not used by \code{sgcbDE})
#' @param adamBeta1 Stored legacy optimizer parameter (currently not used by \code{sgcbDE})
#' @param adamBeta2 Stored legacy optimizer parameter (currently not used by \code{sgcbDE})
#' @param adamEps Stored legacy optimizer parameter (currently not used by \code{sgcbDE})
#' @param weightDecay Stored weight-decay parameter (currently not used by \code{sgcbDE})
#' @param gradClip Stored gradient-clipping threshold (currently not used by \code{sgcbDE})
#' @param lrDecayFactor Stored learning-rate decay factor (currently not used by \code{sgcbDE})
#' @param lrDecayPatience Stored learning-rate decay patience (currently not used by \code{sgcbDE})
#' @param earlyStopPatience Stored early-stopping patience (currently not used by \code{sgcbDE})
#' @param minDelta Stored early-stopping tolerance (currently not used by \code{sgcbDE})
#' @param bootB Number of bootstrap replicates used when \code{bootstrap = TRUE}
#' @param mGrid Stored calibrated-bootstrap grid (currently not used by \code{sgcbDE})
#' @param cbReps Stored calibrated-bootstrap repetitions (currently not used by \code{sgcbDE})
#' @param cbTargetCoverage Stored calibrated-bootstrap target coverage (currently not used by \code{sgcbDE})
#' @param ciLevel Stored confidence-interval level (currently not used by \code{sgcbDE})
#' @param paramLo Stored lower parameter bound (currently not used by \code{sgcbDE})
#' @param paramHi Stored upper parameter bound (currently not used by \code{sgcbDE})
#' @param nDropout Stored dropout replicate count (currently not used by \code{sgcbDE})
#' @param dropoutRate Stored dropout rate (currently not used by \code{sgcbDE})
#' @param nCores Stored parallel-core count (currently not used by \code{sgcbDE})
#' @param seed Stored random seed (currently not used by \code{sgcbDE})
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

  methods::new("SGCBConfig",
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
#' @description Return the default \code{SGCBConfig} object. As with
#' \code{SGCBConfig()}, only the \code{bootB} slot is currently read by the
#' main \code{sgcbDE()} pipeline; the remaining slots are retained for
#' compatibility and future extension.
#' @return SGCBConfig object
#' @export
defaultSGCBConfig <- function() {
  SGCBConfig()
}
