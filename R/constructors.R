# =============================================================================
# SGCB: Constructor functions
# =============================================================================

# -----------------------------------------------------------------------------
# SGCBConfig constructor
# -----------------------------------------------------------------------------
#' @title Create SGCBConfig configuration object
#' @description Create an \code{SGCBConfig} object controlling bootstrap and
#' parallelisation behaviour for \code{sgcbDE()}.
#' @param bootB Number of bootstrap replicates (used when \code{bootstrap = TRUE})
#' @param nCores Number of parallel cores
#' @param seed Random seed for reproducibility
#' @param ... Ignored (backward compatibility with legacy slot names)
#' @return SGCBConfig object
#' @export
SGCBConfig <- function(
  bootB = 1000L,
  nCores = parallel::detectCores(),
  seed = 12345L,
  ...
) {
  methods::new("SGCBConfig",
    bootB = as.integer(bootB),
    nCores = as.integer(nCores),
    seed = as.integer(seed)
  )
}

# -----------------------------------------------------------------------------
# Default configuration
# -----------------------------------------------------------------------------
#' @title Get default configuration
#' @description Return the default \code{SGCBConfig} object.
#' @return SGCBConfig object
#' @export
defaultSGCBConfig <- function() {
  SGCBConfig()
}
