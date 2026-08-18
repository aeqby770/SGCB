# =============================================================================
# 51_extract_null_pvalues.R
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
})

CFG <- list(
  V2_SGCB  = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/sgcb/results"),
  V2_OTHER = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/other_methods/results"),
  OUT      = paste0(PROJECT_ROOT, "/benchmark/output/paper_data_v2/bulk_null_pvalues.csv.gz")
)

OTHER_METHODS <- c("DESeq2", "EBSeq", "NOISeq", "edgeR", "glmGamPoi", "limma", "samr")

collect_pvalues <- function(path, method, replicate) {
  if (!file.exists(path)) {
    message("  missing: ", basename(path))
    return(NULL)
  }
  x <- readRDS(path)
  pv <- x$pvalue
  if (is.null(pv) || length(pv) == 0L) {
    message("  no pvalue column: ", basename(path))
    return(NULL)
  }
  data.table(method = method, replicate = replicate, pvalue = as.numeric(pv))
}

out <- list()
cat("Reading SGCB null_fdr rds (10 files)...\n")
for (i in 1:10) {
  out[[length(out) + 1L]] <- collect_pvalues(
    file.path(CFG$V2_SGCB, sprintf("null_fdr__%d_full.rds", i)),
    "SGCB", i
  )
}
for (m in OTHER_METHODS) {
  cat(sprintf("Reading %s null_fdr rds (10 files)...\n", m))
  for (i in 1:10) {
    out[[length(out) + 1L]] <- collect_pvalues(
      file.path(CFG$V2_OTHER, sprintf("%s__null_fdr__%d_full.rds", m, i)),
      m, i
    )
  }
}

combined <- rbindlist(out, fill = TRUE, use.names = TRUE)
combined <- combined[!is.na(pvalue)]

fwrite(combined, CFG$OUT, compress = "gzip")

cat("\nWrote:", CFG$OUT, "\n")
cat("Total rows:", nrow(combined), "\n")
cat("\nPer-method row counts:\n")
print(combined[, .N, by = method])
