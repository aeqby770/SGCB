# =========================================================================
# 37_v5_run_all.R — SGCB JBHI paper v5 figure regeneration entry point
# 用法: Rscript 37_v5_run_all.R
# 升级要点 (相对 v3):
#   - Fig3/Fig4/Fig6 改用 ComplexHeatmap + row annotation
#   - Fig7 合并原 FigS5 (Cox forest + KM) 解决 reviewer #7
#   - FigS1 manhattan 加强 ggrepel 防重叠 + 饱和色带
#   - FigS2 加 |log2FC|/|log2CV2| 面板注解
#   - FigS4 从 5 panel 精简为 2x2
#   - 输出目录 paper_figures_v5/
# =========================================================================
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
  library(ggtext)
  library(ggrepel)
  library(ggdist)
  library(ggbeeswarm)
  library(ggpointdensity)
  library(ggrastr)
  library(colorspace)
  library(scico)
  library(ragg)
  library(scales)
  library(edgeR)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggsci)
  library(ggpubr)
  library(ggforce)
  library(ggnewscale)
  library(ggExtra)
  library(survival)
  library(survminer)
  library(cowplot)
})

PROJECT_ROOT <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

CONFIG <- list(
  IN_DIR     = file.path(PROJECT_ROOT, "data", "paper_data_v2"),
  REVIEW_DIR = file.path(PROJECT_ROOT, "data", "reviewer_checks"),
  DATA_DIR   = file.path(PROJECT_ROOT, "data"),
  OUT_DIR    = Sys.getenv(
    "SGCB_FIGURE_OUT",
    unset = file.path(PROJECT_ROOT, "results", "regenerated_figures")
  ),
  DPI        = 400L,
  FDR_THRESH = 0.05,
  BASE_SIZE  = 9
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== SGCB JBHI v5 Figure Generation ===\n")
cat("Output:", CONFIG$OUT_DIR, "\n\n")

# --- Load data -----------------------------------------------------------
cat("Loading paper_data_v2 tables...\n")
combined_classic <- fread(file = file.path(CONFIG$IN_DIR, "bulk_classic.csv"))
combined_sim     <- fread(file = file.path(CONFIG$IN_DIR, "bulk_simulation.csv"))
combined_large   <- fread(file = file.path(CONFIG$IN_DIR, "bulk_large_scale.csv"))
combined_null    <- fread(file = file.path(CONFIG$IN_DIR, "bulk_null_fdr_summary.csv"))
combined_null_rep     <- fread(file = file.path(CONFIG$IN_DIR, "bulk_null_fdr.csv"))
combined_null_pvalues <- fread(file = file.path(CONFIG$IN_DIR, "bulk_null_pvalues.csv.gz"))
combined_dv           <- fread(file = file.path(CONFIG$IN_DIR, "bulk_dv_sgcb.csv"))
combined_seqc         <- fread(file = file.path(CONFIG$IN_DIR, "bulk_seqc.csv"))
combined_pb_de        <- fread(file = file.path(CONFIG$IN_DIR, "pseudobulk_de.csv"))
combined_pb_null      <- fread(file = file.path(CONFIG$IN_DIR, "pseudobulk_null.csv"))
combined_prot_de      <- fread(file = file.path(CONFIG$IN_DIR, "proteomics_de.csv"))
combined_prot_null    <- fread(file = file.path(CONFIG$IN_DIR, "proteomics_null.csv"))

cat("Loaded", nrow(combined_classic), "classic,",
    nrow(combined_sim), "sim,",
    nrow(combined_large), "large,",
    nrow(combined_null_pvalues), "null p-values.\n\n")

# --- Run v5 scripts ------------------------------------------------------
cat("=== Part 1: Benchmark (Fig 2-6 + Fig7 SEQC) ===\n")
source(file.path(PROJECT_ROOT, "scripts", "figures", "37b_v5_publication_figures.R"),
        local = FALSE, encoding = "UTF-8")

cat("\n=== Part 2: Biological + Supplementary (Fig 7 TCGA / Fig 8 / FigS1-S5) ===\n")
source(file.path(PROJECT_ROOT, "scripts", "figures", "37c_v5_publication_figures.R"),
        local = FALSE, encoding = "UTF-8")

cat("\n=== Part 3: route-C redesign overlay (overwrites redesigned figures) ===\n")
source(file.path(PROJECT_ROOT, "scripts", "figures", "37d_v5_redesign.R"),
        local = FALSE, encoding = "UTF-8")

cat("\n=== All v5 figures generated ===\n")
cat("Output directory:", CONFIG$OUT_DIR, "\n")
