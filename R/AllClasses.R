# =============================================================================
# SGCB: S4 class definitions
# =============================================================================

# -----------------------------------------------------------------------------
# SGCBConfig: configuration class
# -----------------------------------------------------------------------------
#' @title SGCBConfig class
#' @description Stores SGCB analysis configuration parameters.
#' @slot bootB Number of bootstrap replicates used when \code{bootstrap = TRUE}
#' @slot nCores Number of parallel cores
#' @slot seed Random seed for reproducibility
#' @exportClass SGCBConfig
setClass("SGCBConfig",
  slots = c(
    bootB = "integer",
    nCores = "integer",
    seed = "integer"
  ),
  prototype = list(
    bootB = 1000L,
    nCores = 1L,
    seed = 12345L
  )
)
