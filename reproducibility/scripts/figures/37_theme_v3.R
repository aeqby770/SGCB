
# =========================================================================
# 37_theme_v3.R — SGCB JBHI 论文统一高级可视化主题 + 配色方案
# 所有 37b_v3 / 37c_v3 / story 可视化脚本 source 此文件
# =========================================================================

# --- 1. 配色体系 -----------------------------------------------------------

# 方法配色: 基于 ggsci NPG (Nature Publishing Group) 配色，SGCB 使用标志性红色
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
METHOD_COLORS_V3 <- c(
  SGCB      = "#E64B35",   # NPG Red (主角突出)
  DESeq2    = "#4DBBD5",   # NPG Blue
  edgeR     = "#00A087",   # NPG Green
  limma     = "#3C5488",   # NPG Dark Blue
  NOISeq    = "#F39B7F",   # NPG Orange
  EBSeq     = "#8491B4",   # NPG Purple-ish
  samr      = "#91D1C2",   # NPG Cyan
  glmGamPoi = "#7E6148"    # NPG Brown
)
METHOD_ORDER_V3 <- names(METHOD_COLORS_V3)

# scRNA 全方法配色 (17个, 排除 DiPhiSeq)
SCRNA_COLORS_V3 <- c(
  SGCB            = "#E64B35",
  DESeq2          = "#4DBBD5",
  edgeR           = "#00A087",
  limma           = "#3C5488",
  glmGamPoi       = "#7E6148",
  Seurat          = "#F39B7F",
  MAST            = "#8491B4",
  scDD            = "#91D1C2",
  DEsingle        = "#DC0000",
  distinct        = "#B09C85",
  BPSC            = "#FF7F00",
  muscat          = "#1B9E77",
  aggregateBioVar = "#A6761D",
  glmmTMB         = "#666666",
  NEBULA          = "#E41A1C",
  dv_diffVar      = "#377EB8",
  dv_GAMLSS       = "#984EA3"
)

# Proteomics 全方法配色 (12个, 排除 DiPhiSeq)
PROT_COLORS_V3 <- c(
  SGCB       = "#E64B35",
  limma      = "#3C5488",
  DEqMS      = "#4DBBD5",
  proDA      = "#00A087",
  ROTS       = "#F39B7F",
  msqrob2    = "#8491B4",
  ttest      = "#91D1C2",
  MSstats    = "#7E6148",
  prolfqua   = "#B09C85",
  dv_diffVar = "#377EB8",
  dv_GAMLSS  = "#984EA3"
)

# 合并所有方法配色 (去重)
OMICS_METHOD_COLORS_V3 <- c(METHOD_COLORS_V3, SCRNA_COLORS_V3, PROT_COLORS_V3)
OMICS_METHOD_COLORS_V3 <- OMICS_METHOD_COLORS_V3[!duplicated(names(OMICS_METHOD_COLORS_V3))]

# 方法列表常量
BULK_METHODS_V3 <- c("SGCB", "DESeq2", "edgeR", "limma", "glmGamPoi", "NOISeq", "EBSeq", "samr")
PB_METHODS_V3   <- names(SCRNA_COLORS_V3)
PROT_METHODS_V3 <- names(PROT_COLORS_V3)
DV_METHODS_V3   <- c("SGCB", "dv_diffVar", "dv_GAMLSS", "dv_MDSeq", "dv_clrDV")

# TCGA 通道配色: 更饱和、更有区分度
CHANNEL_COLORS_V3 <- c(
  "DE-only"           = "#3B4CC0",   # 深蓝
  "DV-only"           = "#1B9E77",   # 翡翠绿
  "DE & DV"           = "#D4A017",   # 暗金
  "DD-only (not DE)"  = "#D95F02",   # 橙红
  "DG-only"           = "#7570B3",   # 靛紫
  "overlap"           = "#B3B3B3",   # 中灰
  "Not significant"   = "#E5E5E5"    # 浅灰
)

# FDR 校准专用: 渐变红线
FDR_LINE_COLOR <- "#B22222"

# 渐变色方案
GRADIENT_BLUE_RED <- c("#3B4CC0", "#7B8FD4", "#C4C9E8", "#F7F7F7",
                       "#F0C5AB", "#E07B5E", "#B40426")
GRADIENT_VIRIDIS  <- "viridis"   # 用于 F1 热图
GRADIENT_MAKO     <- "mako"      # 用于 FDR 热图 (scico)

# 生物学故事专用
BIO_GROUP_COLORS <- c(
  normal = "#4DAF4A",
  tumor  = "#E41A1C",
  ctrl   = "#4C72B0",
  stim   = "#C44E52"
)

# --- 2. 主题体系 -----------------------------------------------------------

THEME_JBHI <- function(base_size = 12) {
  # Absolute font sizes derived from base_size (cf. transformGamPoi Nature Methods style)
  fs       <- base_size          # axis title, strip text
  fs_small <- base_size * 0.85   # axis text, legend text
  fs_tiny  <- base_size * 0.75   # caption
  fs_large <- base_size * 1.08   # plot title
  lw       <- 0.3                # universal thin line width

  cowplot::theme_cowplot(
    font_size   = fs,
    rel_small   = fs_small / fs,
    rel_tiny    = fs_tiny / fs,
    rel_large   = fs_large / fs,
    line_size   = lw
  ) %+replace%
    theme(
      # 面板: 细黑框, 无 grid, 白底
      panel.background  = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank(),
      panel.border      = element_rect(color = "grey20", fill = NA, linewidth = lw),
      plot.background   = element_rect(fill = "white", color = NA),
      # 坐标轴: 隐藏线条 (由 panel.border 代替), 细 ticks
      axis.line         = element_blank(),
      axis.text         = element_text(color = "black", size = fs_small),
      axis.title        = element_text(size = fs),
      axis.ticks        = element_line(color = "grey20", linewidth = lw),
      axis.ticks.length = unit(1.2, "mm"),
      # 图例: 紧凑
      legend.position    = "bottom",
      legend.box         = "horizontal",
      legend.key.width   = unit(0.9, "lines"),
      legend.key.height  = unit(0.65, "lines"),
      legend.title       = element_text(face = "bold", size = fs_small),
      legend.text        = element_text(size = fs_tiny),
      legend.background  = element_rect(fill = "white", color = NA),
      legend.key         = element_blank(),
      legend.margin      = margin(t = 1, b = 1),
      legend.spacing.x   = unit(3, "pt"),
      # 分面: 无框白底 (干净), 或极细灰框
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", color = "black",
                                       size = fs_small, margin = margin(t = 2, b = 2)),
      panel.spacing    = unit(0.4, "lines"),
      # 标题: ggtext markdown 支持
      plot.title         = ggtext::element_markdown(face = "bold", size = fs_large,
                                                     lineheight = 1.15, margin = margin(b = 4),
                                                     hjust = 0),
      plot.subtitle      = ggtext::element_markdown(size = fs_small, color = "#444444",
                                                     margin = margin(b = 4), hjust = 0),
      plot.caption       = element_text(size = fs_tiny, color = "#666666", hjust = 0),
      plot.title.position = "plot",
      plot.margin        = margin(6, 8, 4, 8)
    )
}

# 热图专用: 去掉边框, 只保留 tile
THEME_HEATMAP <- function(base_size = 11) {
  THEME_JBHI(base_size) %+replace%
    theme(
      panel.border = element_blank(),
      axis.ticks   = element_blank(),
      axis.line    = element_blank()
    )
}

# --- ggplot2 geom defaults: 与 theme 线宽一致 (cf. transformGamPoi) ---
update_geom_defaults("line",    list(linewidth = 0.4))
update_geom_defaults("hline",   list(linewidth = 0.3))
update_geom_defaults("vline",   list(linewidth = 0.3))
update_geom_defaults("segment", list(linewidth = 0.3))
update_geom_defaults("col",     list(linewidth = 0.2))
update_geom_defaults("bar",     list(linewidth = 0.2))
update_geom_defaults("point",   list(size = 1.5))
update_geom_defaults("text",    list(size = 3))

# --- 3. 保存函数 -----------------------------------------------------------

SAVE_FIG_V3 <- function(nm, pl, w, h, out_dir = CONFIG$OUT_DIR, dpi = 400L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(nm, ".pdf")), pl,
         width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(nm, ".png")), pl,
         width = w, height = h, dpi = dpi, device = ragg::agg_png, bg = "white")
  cat("  Saved:", nm, "(", w, "x", h, ")\n")
}

# --- 4. 工具函数 -----------------------------------------------------------

# 自动 best-in-bold: 给 data.table 每列的最大值加粗 markdown
bold_best <- function(dt, value_cols, direction = "max") {
  dt2 <- copy(dt)
  for (vc in value_cols) {
    vals <- dt2[[vc]]
    if (!is.numeric(vals)) next
    best_idx <- if (direction == "max") which.max(vals) else which.min(vals)
    formatted <- sprintf("%.3f", vals)
    formatted[best_idx] <- paste0("**", formatted[best_idx], "**")
    dt2[, (vc) := formatted]
  }
  dt2
}

# 分面标签清理
clean_dataset_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  x
}

cat("[theme_v3] Loaded SGCB JBHI publication theme.\n")
