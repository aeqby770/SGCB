# =========================================================================
# 37_theme_v5.R — SGCB JBHI paper v5 可视化主题 + PLOTR 集成辅助函数
# 被 37b_v5 / 37c_v5 source 进来
# 设计要点:
#   - 沿用 v3 的 METHOD_COLORS_V3 / CHANNEL_COLORS_V3 / THEME_JBHI() 命名
#   - 用 ggsci NPG 作为主色，但对 channel 另起一套更深的学术配色
#   - 提供 ComplexHeatmap 封装: ch_simple_heatmap() 和 ch_with_row_ann()
#   - 提供 ggrastr 封装和公用保存函数 SAVE_FIG_V5
# =========================================================================
source(file.path(PROJECT_ROOT, "scripts", "figures", "37_theme_v3.R"),
       local = FALSE, encoding = "UTF-8")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggsci)
  library(ggpubr)
  library(ggforce)
  library(ggrastr)
  library(ggrepel)
  library(ggtext)
  library(ggnewscale)
  library(ggdist)
  library(patchwork)
  library(cowplot)
  library(colorspace)
  library(scales)
})

# --- v5 high-impact palette override ---------------------------------------
# Keep the v3 object names so existing scripts stay compatible, but move away
# from the default NPG look. The palette is intentionally higher contrast for
# bioinformatics method figures where first-impression signal matters.
METHOD_COLORS_V3 <- c(
  SGCB      = "#FF3B30",
  DESeq2    = "#00B7E5",
  edgeR     = "#00B894",
  limma     = "#3557B7",
  NOISeq    = "#FF8A65",
  EBSeq     = "#8A9CC8",
  samr      = "#76D7C4",
  glmGamPoi = "#8A6F52"
)
METHOD_ORDER_V3 <- names(METHOD_COLORS_V3)

SCRNA_COLORS_V3 <- c(
  SGCB            = "#FF3B30",
  DESeq2          = "#00B7E5",
  edgeR           = "#00B894",
  limma           = "#3557B7",
  glmGamPoi       = "#8A6F52",
  Seurat          = "#BDBDBD",
  MAST            = "#6F86C6",
  scDD            = "#39C7B8",
  DEsingle        = "#F9704F",
  distinct        = "#A58D78",
  BPSC            = "#FFB000",
  muscat          = "#23A56F",
  aggregateBioVar = "#B37A22",
  glmmTMB         = "#5C5C68",
  NEBULA          = "#E0316B",
  dv_diffVar      = "#2A6FDB",
  dv_GAMLSS       = "#A64AC9"
)

PROT_COLORS_V3 <- c(
  SGCB       = "#FF3B30",
  limma      = "#3557B7",
  DEqMS      = "#00B7E5",
  proDA      = "#00B894",
  ROTS       = "#FF8A65",
  msqrob2    = "#8A9CC8",
  ttest      = "#76D7C4",
  MSstats    = "#8A6F52",
  prolfqua   = "#B8A28A",
  dv_diffVar = "#2A6FDB",
  dv_GAMLSS  = "#A64AC9"
)

OMICS_METHOD_COLORS_V3 <- c(METHOD_COLORS_V3, SCRNA_COLORS_V3, PROT_COLORS_V3)
OMICS_METHOD_COLORS_V3 <- OMICS_METHOD_COLORS_V3[!duplicated(names(OMICS_METHOD_COLORS_V3))]

CHANNEL_COLORS_V3 <- c(
  "DE-only"           = "#3547C8",
  "DV-only"           = "#00A77E",
  "DG-only"           = "#7E57C2",
  "DE & DV"           = "#D9A20B",
  "DD-only (not DE)"  = "#FF6D00",
  "overlap"           = "#7A7A86",
  "Not significant"   = "#E6E8EF"
)

# --- 继承 v3 配色，但对热图用更正式的色阶 ---
HEATMAP_F1_COLS  <- colorRamp2(c(0, 0.5, 1),
                                c("#F7F9FC", "#89B3D9", "#1F4E79"))
HEATMAP_FDR_COLS <- colorRamp2(c(0, 0.05, 0.15, 0.5),
                                c("#F7F9FC", "#FFE1C7", "#E07B5E", "#8B0000"))
HEATMAP_COUNT_COLS <- colorRamp2(c(0, 1, 2.5, 4.5),
                                  c("#F7F9FC", "#C4DFE6", "#4C72B0", "#1A3A61"))
HEATMAP_CORR_COLS <- colorRamp2(c(-1, 0, 1),
                                 c("#3B4CC0", "white", "#B40426"))

# 面向 Nature/JBHI 的字号统一 (pt) — 与 THEME_JBHI(base_size=10) 对齐
# base=10 → fs_small=8.5, fs_tiny=7.5, fs_large=10.8
FS <- list(title=10.8, sub=8.5, axis=8.5, tick=7.5, legend=7.5, cell=6.5)

# 标准化 colorbar 尺寸 (避免散落各处不一致)
CBAR_H <- guide_colorbar(barwidth = 0.4, barheight = 5, title.position = "top")
CBAR_W <- guide_colorbar(barwidth = 8,   barheight = 0.35, title.position = "left",
                          title.vjust = 1)

# --- ComplexHeatmap 封装：带 legend、grid lines、cell text ---
ch_simple <- function(mat, col_fun, name = "",
                       cell_fmt = "%.2f",
                       row_title = NULL, column_title = NULL,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       row_names_side = "left",
                       column_names_side = "bottom",
                       column_names_rot = 45,
                       rect_gp = gpar(col = "white", lwd = 1),
                       row_title_rot = 0,
                       show_heatmap_legend = TRUE,
                       legend_title_position = "topleft",
                       ...) {
  ht <- ComplexHeatmap::Heatmap(
    mat, col = col_fun, name = name,
    cluster_rows = cluster_rows, cluster_columns = cluster_cols,
    row_names_side = row_names_side,
    column_names_side = column_names_side,
    column_names_rot = column_names_rot,
    rect_gp = rect_gp,
    row_title = row_title, column_title = column_title,
    row_title_rot = row_title_rot,
    row_names_gp = gpar(fontsize = FS$axis),
    column_names_gp = gpar(fontsize = FS$axis),
    cell_fun = if (is.null(cell_fmt)) NULL else function(j, i, x, y, w, h, fill) {
      val <- mat[i, j]
      if (is.na(val)) return(invisible())
      txt <- sprintf(cell_fmt, val)
      # 自动选字色
      lum <- colSums(col2rgb(fill) * c(0.299, 0.587, 0.114))
      col <- if (lum < 140) "white" else "black"
      grid.text(txt, x, y, gp = gpar(fontsize = FS$cell, col = col))
    },
    show_heatmap_legend = show_heatmap_legend,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = FS$legend, fontface = "bold"),
      labels_gp = gpar(fontsize = FS$legend),
      grid_width = unit(3.5, "mm"),
      title_position = legend_title_position
    ),
    ...
  )
  ht
}

# --- 保存函数 ---
SAVE_FIG_V5 <- function(nm, pl, w, h, out_dir = CONFIG$OUT_DIR,
                         dpi = 400L, fmt = c("png", "pdf")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if ("pdf" %in% fmt) {
    ggplot2::ggsave(file.path(out_dir, paste0(nm, ".pdf")), pl,
                     width = w, height = h, device = cairo_pdf, bg = "white")
  }
  if ("png" %in% fmt) {
    ggplot2::ggsave(file.path(out_dir, paste0(nm, ".png")), pl,
                     width = w, height = h, dpi = dpi,
                     device = ragg::agg_png, bg = "white")
  }
  cat("  Saved:", nm, "(", w, "x", h, ")\n")
}

# 保存 ComplexHeatmap (grid graphics) —— 接受 grob 或 HeatmapList
SAVE_HT_V5 <- function(nm, ht, w, h, out_dir = CONFIG$OUT_DIR,
                        dpi = 400L, merge_legend = TRUE,
                        heatmap_legend_side = "right",
                        annotation_legend_side = "right") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  render <- function() {
    if (inherits(ht, c("Heatmap", "HeatmapList"))) {
      ComplexHeatmap::draw(ht,
        merge_legend = merge_legend,
        heatmap_legend_side = heatmap_legend_side,
        annotation_legend_side = annotation_legend_side,
        padding = unit(c(4, 4, 4, 4), "mm"))
    } else if (inherits(ht, "gTree") || inherits(ht, "grob")) {
      grid::grid.draw(ht)
    } else {
      stop("SAVE_HT_V5: unsupported object class ", class(ht)[1])
    }
  }
  png_path <- file.path(out_dir, paste0(nm, ".png"))
  ragg::agg_png(png_path, width = w, height = h, units = "in",
                 res = dpi, bg = "white")
  render(); dev.off()
  pdf_path <- file.path(out_dir, paste0(nm, ".pdf"))
  cairo_pdf(pdf_path, width = w, height = h, bg = "white")
  render(); dev.off()
  cat("  Saved(HT):", nm, "(", w, "x", h, ")\n")
}

# --- grid -> ggplot 转换 (用于 patchwork 拼图) ---
ht_to_grob <- function(ht) {
  grid::grid.grabExpr(ComplexHeatmap::draw(ht,
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"))
}

# --- 浅色标题 span（替代重复的 plot_annotation 手工样式） ---
fig_title <- function(title, subtitle = NULL) {
  patchwork::plot_annotation(
    title = title,
    subtitle = subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = FS$title,
                                  margin = margin(b = 3)),
      plot.subtitle = element_text(size = FS$sub, color = "grey35",
                                     margin = margin(b = 4))))
}

cat("[theme_v5] v5 theme + PLOTR helpers loaded.\n")
