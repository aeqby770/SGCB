# =============================================================================
# 52_kirc_channels_and_go.R
#
#   - paper_data_v2/tcga_kirc_sgcb_channels.csv
#   - paper_data_v2/tcga_pan_cancer_sgcb_channels.csv
#   - paper_data_v2/tcga_kirc_go_de_only.csv
#   - paper_data_v2/tcga_kirc_go_dv_only.csv
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

CFG <- list(
  V2_SGCB = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/sgcb/results"),
  OUT_DIR = paste0(PROJECT_ROOT, "/benchmark/output/paper_data_v2"),
  FDR     = 0.05,
  CANCERS = c("brca", "coad", "hnsc", "kirc", "lihc", "luad", "stad"),
  FOCAL   = "kirc",
  GO_Q    = 0.20,
  GO_P    = 0.05
)

v2_to_legacy <- function(x) {
  dt <- data.table(
    gene                = as.character(x$gene_id),
    padj_de             = as.numeric(x$padj),
    padj_dv             = as.numeric(x$dv_padj),
    padj_dg_a           = as.numeric(x$dg_alpha_padj),
    padj_dg_g           = as.numeric(x$dg_gamma_padj),
    padj_dd             = as.numeric(x$padj_dd),
    log2FC              = as.numeric(x$log2FoldChange),
    dv_log2_cv2_ratio   = as.numeric(x$dv_log2_cv2_ratio),
    dg_gamma_log2_ratio = as.numeric(x$dg_gamma_log2_ratio),
    sgcb_score          = as.numeric(x$SGCB_Score)
  )
  fdr <- CFG$FDR
  dt[, sig_de := !is.na(padj_de) & padj_de < fdr]
  dt[, sig_dv := !is.na(padj_dv) & padj_dv < fdr]
  dt[, sig_dd := !is.na(padj_dd) & padj_dd < fdr]
  dt[, cat_de_only := sig_de & !sig_dv]
  dt[, cat_dv_only := !sig_de &  sig_dv]
  dt[, cat_both    :=  sig_de &  sig_dv]
  dt[, cat_dd_only := !sig_de &  sig_dd]
  dt
}

cat("== Pan-cancer channel summary ==\n")
pan_rows <- list()
focal_dt <- NULL
for (cc in CFG$CANCERS) {
  rds_path <- file.path(CFG$V2_SGCB, sprintf("tcga__%s_full.rds", cc))
  if (!file.exists(rds_path)) {
    message("  missing: ", rds_path); next
  }
  raw  <- readRDS(rds_path)
  dtcc <- v2_to_legacy(raw)
  pan_rows[[cc]] <- data.table(
    cancer     = toupper(cc),
    n_genes    = nrow(dtcc),
    de_only    = sum(dtcc$cat_de_only),
    dv_only    = sum(dtcc$cat_dv_only),
    both_de_dv = sum(dtcc$cat_both),
    dd_only    = sum(dtcc$cat_dd_only)
  )
  if (cc == CFG$FOCAL) focal_dt <- dtcc
  cat(sprintf("  %s: n=%d  DE-only=%d  DV-only=%d  DE&DV=%d  DD-only=%d\n",
              toupper(cc), nrow(dtcc),
              sum(dtcc$cat_de_only), sum(dtcc$cat_dv_only),
              sum(dtcc$cat_both),    sum(dtcc$cat_dd_only)))
}

pan_dt <- rbindlist(pan_rows, use.names = TRUE)
fwrite(pan_dt, file.path(CFG$OUT_DIR, "tcga_pan_cancer_sgcb_channels.csv"))
cat("Pan-cancer CSV ->", file.path(CFG$OUT_DIR, "tcga_pan_cancer_sgcb_channels.csv"), "\n")

# Focal (KIRC) per-gene channels CSV
stopifnot(!is.null(focal_dt))
fwrite(focal_dt, file.path(CFG$OUT_DIR, "tcga_kirc_sgcb_channels.csv"))
cat("KIRC per-gene CSV ->", file.path(CFG$OUT_DIR, "tcga_kirc_sgcb_channels.csv"),
    "(", nrow(focal_dt), "genes )\n")

# --- Step 3: KIRC GO:BP enrichment for DE-only and DV-only ---
cat("\n== KIRC GO:BP enrichment ==\n")

de_only_genes <- focal_dt[cat_de_only == TRUE, gene]
dv_only_genes <- focal_dt[cat_dv_only == TRUE, gene]
bg_genes      <- focal_dt[, gene]
cat(sprintf("  DE-only foreground: %d  DV-only foreground: %d  background: %d\n",
            length(de_only_genes), length(dv_only_genes), length(bg_genes)))

run_go <- function(fg, bg, label) {
  if (length(fg) < 10L) {
    message("  skip ", label, " (fg too small: ", length(fg), ")")
    return(NULL)
  }
  res <- clusterProfiler::enrichGO(
    gene          = fg,
    universe      = bg,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = CFG$GO_P,
    qvalueCutoff  = CFG$GO_Q,
    readable      = FALSE
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0L) {
    message("  ", label, ": 0 enriched terms")
    return(data.table(ID = character(0), Description = character(0),
                      GeneRatio = character(0), BgRatio = character(0),
                      pvalue = numeric(0), p.adjust = numeric(0),
                      qvalue = numeric(0), geneID = character(0),
                      Count = integer(0)))
  }
  dt <- as.data.table(as.data.frame(res))
  cat(sprintf("  %s: %d enriched GO:BP terms (q<%s, p<%s)\n",
              label, nrow(dt), CFG$GO_Q, CFG$GO_P))
  dt
}

go_de <- run_go(de_only_genes, bg_genes, "DE-only")
go_dv <- run_go(dv_only_genes, bg_genes, "DV-only")

if (!is.null(go_de))
  fwrite(go_de, file.path(CFG$OUT_DIR, "tcga_kirc_go_de_only.csv"))
if (!is.null(go_dv))
  fwrite(go_dv, file.path(CFG$OUT_DIR, "tcga_kirc_go_dv_only.csv"))

cat("\nGO CSV ->\n  ",
    file.path(CFG$OUT_DIR, "tcga_kirc_go_de_only.csv"), "\n  ",
    file.path(CFG$OUT_DIR, "tcga_kirc_go_dv_only.csv"), "\n")

cat("\n== done ==\n")
