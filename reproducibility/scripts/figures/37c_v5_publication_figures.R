# =========================================================================
# 37c_v5 — SGCB JBHI paper Fig 7 (TCGA+survival association merged) / Fig 8 (kang18)
#          + FigS1 (manhattan) / FigS2 (gene boxplots) / FigS3 (pan-cancer)
#          / FigS4 (reviewer) / FigS5 (clinical standalone, backup)
# 数据由 37_v5_run_all.R 入口加载 (combined_*)
# 关键升级:
#   - Fig7 新增 Cox forest (C) 和 KM (D) 两个 survival-association 面板
#   - FigS1 manhattan 用 ggrepel + ggforce::geom_mark_rect 防重叠
#   - FigS2 每格加 |log2FC| / |log2CV2| 双注解
#   - FigS3 stacked bar 标签阈值降低
#   - FigS4 从 5 panel 精简为 2x2 布局
# =========================================================================

source(file.path(PROJECT_ROOT, "scripts", "figures", "37_theme_v5.R"),
       local = FALSE, encoding = "UTF-8")
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggExtra)
})

# -- 共享资源 --------------------------------------------------------------
kirc_dt_all <- fread(file = file.path(CONFIG$IN_DIR, "tcga_kirc_sgcb_channels.csv"))
kirc_dt_all[, sig_dg := pmin(padj_dg_a, padj_dg_g, na.rm = TRUE) < 0.05]
kirc_dt_all[, cat_dg_only := sig_dg & !sig_de & !sig_dv]
kirc_dt_all[, channel := fifelse(cat_both == TRUE, "DE & DV",
  fifelse(cat_de_only == TRUE, "DE-only",
  fifelse(cat_dv_only == TRUE, "DV-only",
  fifelse(cat_dg_only == TRUE, "DG-only",
  fifelse(cat_dd_only == TRUE, "DD-only (not DE)", "Not significant")))))]
CH_LEVELS <- c("DE-only", "DV-only", "DG-only", "DE & DV",
                "DD-only (not DE)", "Not significant")
kirc_dt_all[, channel := factor(channel, levels = CH_LEVELS)]

CH_COLS_V5 <- c(
  "DE-only"          = "#2F5BFF",
  "DV-only"          = "#15C8A2",
  "DG-only"          = "#8A4DFF",
  "DE & DV"          = "#FFB000",
  "DD-only (not DE)" = "#FF5A1F",
  "Not significant"  = "#D0D4DB")

# =========================================================================
# Run 51_tcga_kirc_clinical_extended.R to obtain cox_all / matched / km_fit
# =========================================================================
cat("  Sourcing 51_tcga_kirc_clinical_extended.R for survival data...\n")
local({
  env <- new.env()
  env$FIG_DIR_OVERRIDE <- CONFIG$OUT_DIR
  source(file.path(PROJECT_ROOT, "scripts", "figures", "51_tcga_kirc_clinical_extended.R"),
          local = env, encoding = "UTF-8")
  assign("SURV_ENV", env, envir = .GlobalEnv)
})
matched    <- SURV_ENV$matched
cox_all    <- SURV_ENV$cox_all
stage_tab  <- SURV_ENV$stage_tab
gene_sets  <- SURV_ENV$gene_sets
channel_label <- SURV_ENV$channel_label

# =========================================================================
# Fig 7 — TCGA-KIRC merged (channel comp + landscape + Cox + KM + violins)
# =========================================================================
cat("  Fig 7 (v5): TCGA-KIRC merged (4 panels)...\n")

# --- Panel A: channel composition ---------------------------------------
active_total <- sum(kirc_dt_all$channel != "Not significant")
ch_counts <- kirc_dt_all[channel != "Not significant",
                           .N, by = channel][order(channel)]
ch_counts[, frac := N / active_total]
ch_counts[, end := cumsum(frac) * 2 * pi]
ch_counts[, start := shift(end, fill = 0)]
ch_counts[, mid := (start + end) / 2]
ch_counts[, x_lab := 0.74 * cos(mid)]
ch_counts[, y_lab := 0.74 * sin(mid)]
ch_counts[, x_callout := 1.05 * cos(mid)]
ch_counts[, y_callout := 1.05 * sin(mid)]
ch_counts[, channel_short := fifelse(as.character(channel) == "DD-only (not DE)",
                                      "DD-only", as.character(channel))]
ch_counts[, label := sprintf("%s\n%s", channel_short, scales::comma(N))]
p7a <- ggplot(ch_counts) +
  ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.44, r = 1.0,
                            start = start, end = end, fill = channel),
                        color = "white", linewidth = 0.72,
                        show.legend = FALSE) +
  geom_text(data = ch_counts[frac >= 0.02],
            aes(x = x_lab, y = y_lab, label = label),
            size = 2.35, lineheight = 0.82, fontface = "bold",
            color = "#15151A") +
  ggrepel::geom_label_repel(data = ch_counts[frac < 0.02],
            aes(x = x_callout, y = y_callout, label = label, fill = channel),
            color = "#15151A", size = 2.1, fontface = "bold",
            label.size = 0.18, label.padding = unit(0.07, "lines"),
            min.segment.length = 0, segment.color = "#606876",
            segment.size = 0.18, seed = 42, show.legend = FALSE) +
  annotate("text", x = 0, y = 0.08, label = scales::comma(active_total),
           size = 4.8, fontface = "bold", color = "#111827") +
  annotate("text", x = 0, y = -0.12, label = "active\ngenes",
           size = 2.35, lineheight = 0.78, color = "#4B5563") +
  scale_fill_manual(values = CH_COLS_V5) +
  coord_fixed(xlim = c(-1.12, 1.12), ylim = c(-1.12, 1.12), clip = "off") +
  labs(title = "**A.**",
        x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(plot.title = ggtext::element_markdown(size = 11, face = "bold",
                                               hjust = 0, lineheight = 1.0),
        plot.margin = margin(3, 8, 3, 8))

# --- Panel B: channel landscape scatter with marginal densities ----------
kirc_act <- kirc_dt_all[channel != "Not significant"]
kirc_act[, abs_lfc := pmin(abs(log2FC), quantile(abs(log2FC), 0.995,
                                                     na.rm = TRUE))]
kirc_act[, abs_cv2 := pmin(abs(dv_log2_cv2_ratio),
                              quantile(abs(dv_log2_cv2_ratio), 0.995,
                                         na.rm = TRUE))]
# top 1 per channel by composite score for labels
kirc_act[, score := abs_lfc * abs_cv2]
lbl_genes <- kirc_act[, .SD[order(-score)][1L:2L], by = channel][!is.na(gene)]
lfc_ref <- quantile(kirc_act$abs_lfc, 0.75, na.rm = TRUE)
cv2_ref <- quantile(kirc_act$abs_cv2, 0.75, na.rm = TRUE)

p7b_base <- ggplot(kirc_act, aes(x = abs_lfc, y = abs_cv2)) +
  stat_density_2d(color = "#111827", linewidth = 0.18, alpha = 0.22,
                  bins = 8, show.legend = FALSE) +
  geom_vline(xintercept = lfc_ref, color = "#A3AAB8", linewidth = 0.25,
             linetype = "dashed") +
  geom_hline(yintercept = cv2_ref, color = "#A3AAB8", linewidth = 0.25,
             linetype = "dashed") +
  rasterise(geom_point(aes(color = channel), size = 0.32, alpha = 0.42,
                       shape = 16), dpi = 300) +
  geom_point(data = lbl_genes, aes(fill = channel), shape = 21,
              color = "black", size = 2.45, stroke = 0.55) +
  ggrepel::geom_text_repel(data = lbl_genes, aes(label = gene),
              size = 2.5, fontface = "bold", color = "grey10",
              bg.color = "white", bg.r = 0.15,
              box.padding = 0.3, force = 5, max.overlaps = Inf,
              segment.size = 0.2, segment.color = "grey40", seed = 42,
              show.legend = FALSE) +
  scale_color_manual(values = CH_COLS_V5,
                      breaks = setdiff(CH_LEVELS, "Not significant"),
                      name = NULL) +
  scale_fill_manual(values = CH_COLS_V5, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "**B.**",
        x = expression("| " ~ log[2] ~ "FC |"),
        y = expression("| " ~ log[2] ~ "CV"^2 ~ "ratio |")) +
  THEME_JBHI(base_size = 10) +
  theme(legend.position = "none",
         legend.key.size = unit(0.3, "cm"),
         legend.text = element_text(size = 7),
         panel.background = element_rect(fill = "#FAFBFD", color = NA),
         panel.border = element_rect(color = "#15151A", fill = NA,
                                     linewidth = 0.35))
# Since ggMarginal would collapse legend, skip it; use facet-less density as
# insets would be added by patchwork instead.
p7b <- p7b_base

# --- Panel C: Cox forest plot (from sourced cox_all) ---------------------
cox_plot <- copy(cox_all)
cox_plot[, signature := factor(signature, levels = rev(unname(channel_label)))]
cox_plot[, adjustment := factor(adjustment,
  levels = c("Univariate", "Adjusted (stage + grade + age)"))]
cox_plot[, adjustment_short := factor(
  fifelse(adjustment == "Univariate", "Univ.", "Adjusted"),
  levels = c("Univ.", "Adjusted"))]
cox_plot[, lab := sprintf("HR %.2f\np=%s",
                            HR, formatC(p, format = "g", digits = 2))]
x_hi <- max(cox_plot$upper, na.rm = TRUE) * 2.5
label_x <- max(cox_plot$upper, na.rm = TRUE) * 1.1

p7c <- ggplot(cox_plot, aes(x = HR, y = signature, color = adjustment_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55",
              linewidth = 0.35) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                   position = position_dodge(width = 0.58),
                   size = 0.4, linewidth = 0.5) +
  geom_text(aes(x = label_x, label = lab),
             position = position_dodge(width = 0.58),
             hjust = 0, size = 2.0, lineheight = 0.85,
             show.legend = FALSE) +
  scale_color_manual(values = c("Univ." = "#3B4CC0",
                                 "Adjusted" = "#D4A017"),
                      name = NULL) +
  scale_x_continuous(trans = "log2",
                      limits = c(0.5, x_hi),
                      breaks = c(0.5, 1, 1.5, 2, 3)) +
  labs(title = "**C.**",
        x = expression("Hazard ratio (" ~ log[2] ~ " scale)"),
        y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(legend.position = "bottom",
         legend.key.size = unit(0.3, "cm"),
         legend.text = element_text(size = 7),
         panel.grid.major.y = element_blank(),
         axis.text.y = element_text(face = "bold", size = 8))

# --- Panel D: KM curve for DV-only signature -----------------------------
km_fit <- survfit(Surv(os_time, os_event) ~ group_DV_only, data = matched)
km_obj <- ggsurvplot(km_fit, data = matched,
                       palette = c("#C44E52", "#3B4CC0"),
                       legend = "top",
                       legend.title = "DV-only",
                       legend.labs = c("High", "Low"),
                       pval = TRUE, pval.size = 2.8,
                       pval.coord = c(100, 0.08),
                       risk.table = FALSE, conf.int = TRUE,
                       conf.int.alpha = 0.10,
                       censor.size = 2, size = 0.6,
                       xlab = "Time (days)", ylab = "Overall survival",
                       ggtheme = THEME_JBHI(base_size = 10) +
                         theme(legend.key.size = unit(0.3, "cm"),
                                legend.text = element_text(size = 7),
                                legend.title = element_text(size = 7.5,
                                                              face = "bold")))
p7d <- km_obj$plot +
  labs(title = "**D.**")

# --- Panel E: representative gene violins (6 genes) ---------------------
kirc_expr <- readRDS(file.path(PROJECT_ROOT, "data", "tcga_kirc_cache",
                               "tcga_kirc_tumor_vs_normal.rds"))
counts_kirc <- kirc_expr$counts
group_kirc  <- kirc_expr$group

de_cand <- kirc_dt_all[cat_de_only == TRUE & padj_de < 1e-5 & padj_dv > 0.1][order(-abs(log2FC))]
dv_cand <- kirc_dt_all[cat_dv_only == TRUE & padj_de > 0.3 &
                         abs(log2FC) < 0.5][order(padj_dv)]
dg_cand <- kirc_dt_all[cat_dg_only == TRUE &
                         abs(log2FC) < 0.6 & padj_de > 0.3][order(-abs(dg_gamma_log2_ratio))]
rep_genes <- unique(c(head(de_cand, 2)$gene, head(dv_cand, 2)$gene,
                       head(dg_cand, 2)$gene))
if (length(rep_genes) < 6L) {
  # fallback pool: top ranked per channel by padj
  rep_genes <- unique(c(rep_genes,
    head(kirc_dt_all[cat_de_only == TRUE][order(padj_de)], 2)$gene,
    head(kirc_dt_all[cat_dv_only == TRUE][order(padj_dv)], 2)$gene))
  rep_genes <- rep_genes[!is.na(rep_genes)][1:6]
}

library(edgeR)
dge_b <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts_kirc,
                                                  group = group_kirc))
sf_b <- dge_b$samples$lib.size * dge_b$samples$norm.factors
sf_b <- sf_b / exp(mean(log(sf_b)))
nc_b <- t(t(counts_kirc) / sf_b) + 0.5
expr_mat <- log2(nc_b[rep_genes, , drop = FALSE])

gene_meta <- kirc_dt_all[gene %in% rep_genes]
gene_meta[, channel_simple := fifelse(cat_de_only == TRUE, "DE-only",
  fifelse(cat_dv_only == TRUE, "DV-only",
  fifelse(cat_dg_only == TRUE, "DG-only", "Other")))]
gene_meta[, facet_lbl := paste0(gene, " (", channel_simple, ")")]
gene_meta[, lbl := sprintf("|log2FC|=%.2f, |\u0394log2CV\u00b2|=%.2f",
                              abs(log2FC), abs(dv_log2_cv2_ratio))]
gene_meta[, gene_order := match(gene, rep_genes)]

expr_long <- data.table(
  gene = rep(rep_genes, each = ncol(expr_mat)),
  expr = as.vector(t(expr_mat)),
  group = rep(as.character(group_kirc), length(rep_genes)))
expr_long <- merge(expr_long, gene_meta[, .(gene, facet_lbl, lbl,
                                              channel_simple, gene_order)],
                    by = "gene")
expr_long[, facet_lbl := reorder(facet_lbl, gene_order)]

p7e <- ggplot(expr_long, aes(x = group, y = expr, fill = group)) +
  ggdist::stat_halfeye(adjust = 0.8, width = 0.5, .width = 0,
                         justification = -0.15, point_colour = NA,
                         alpha = 0.5) +
  geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 0.9,
                color = "grey25", linewidth = 0.3) +
  facet_wrap(~facet_lbl, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = c(normal = "#4DAF4A", tumor = "#E41A1C"),
                     guide = "none") +
  geom_text(data = unique(gene_meta[, .(facet_lbl, lbl)]),
             aes(x = 1.5, y = Inf, label = lbl),
             inherit.aes = FALSE, vjust = 1.4, size = 2.0,
             color = "grey25", fontface = "italic") +
  labs(title = "**E.**",
        x = NULL, y = expression(log[2] ~ "expression")) +
  THEME_JBHI(base_size = 10) +
  theme(strip.text = element_text(size = 7.5, face = "bold"),
         axis.text.x = element_text(size = 7),
         panel.spacing = unit(0.3, "lines"))

# --- Fig 7 assembly ------------------------------------------------------
row1 <- cowplot::plot_grid(p7a, p7b, p7c, ncol = 3,
                             rel_widths = c(1, 1.2, 1.3),
                             align = "h", axis = "tb")
row2 <- cowplot::plot_grid(p7d, ncol = 1, align = "h", axis = "tb")
fig7 <- cowplot::plot_grid(row1, row2,
                             ncol = 1, rel_heights = c(1, 0.9))
SAVE_FIG_V5("Fig7_tcga_kirc_merged", fig7, 7.5, 7.9)

# Back-up standalone FigS5 (Cox forest + KM + stage) re-render to v5 dir
if (!is.null(SURV_ENV$fig_combined)) {
  SAVE_FIG_V5("FigS5_tcga_kirc_clinical", SURV_ENV$fig_combined, 12, 8)
}

# =========================================================================
# Fig 8 — PB_kang18 merged (4 panels)
# =========================================================================
cat("  Fig 8 (v5): PB_kang18 merged...\n")
KANG_TBL <- file.path(PROJECT_ROOT, "data", "biological_story_outputs",
                      "pb_kang18", "tables")
kang_ch <- fread(file = file.path(KANG_TBL, "pb_kang18_channel_summary.csv"))
kang_en <- fread(file = file.path(KANG_TBL, "pb_kang18_enrichment.csv"))
kang_mk <- fread(file = file.path(KANG_TBL, "pb_kang18_celltype_marker_enrichment.csv"))
kang_dist <- fread(file = file.path(KANG_TBL, "pb_kang18_representative_distributions_long.csv"))

# --- Panel A: channel counts --------------------------------------------
kang_ch_plot <- kang_ch[channel_category %in% c("DE-only", "DV-only",
                                                    "DG-only", "overlap")]
kang_ch_plot <- kang_ch_plot[, .(n_features = sum(n_features, na.rm = TRUE)),
                              by = channel_category]
kang_ch_plot[, channel_category := factor(channel_category,
  levels = c("DE-only", "DV-only", "DG-only", "overlap"))]
total_feat <- sum(kang_ch_plot$n_features)
kang_ch_plot[, frac := n_features / total_feat]
kang_ch_plot[, end := cumsum(frac) * 2 * pi]
kang_ch_plot[, start := shift(end, fill = 0)]
kang_ch_plot[, mid := (start + end) / 2]
kang_ch_plot[, x_lab := 0.72 * cos(mid)]
kang_ch_plot[, y_lab := 0.72 * sin(mid)]
kang_ch_plot[, x_callout := 1.04 * cos(mid)]
kang_ch_plot[, y_callout := 1.04 * sin(mid)]
kang_ch_plot[, channel_short := fifelse(as.character(channel_category) == "overlap",
                                         "DE & DV", as.character(channel_category))]
kang_ch_plot[, lbl := sprintf("%s\n%s", channel_short, scales::comma(n_features))]
KANG_CH_COLS <- c(CH_COLS_V5[c("DE-only", "DV-only", "DG-only")],
                  overlap = unname(CH_COLS_V5["DE & DV"]))
p8a <- ggplot(kang_ch_plot) +
  ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.42, r = 1.0,
                            start = start, end = end, fill = channel_category),
                        color = "white", linewidth = 0.7,
                        show.legend = FALSE) +
  geom_text(data = kang_ch_plot[frac >= 0.035],
            aes(x = x_lab, y = y_lab, label = lbl),
            size = 2.35, lineheight = 0.82, fontface = "bold",
            color = "#15151A") +
  ggrepel::geom_label_repel(data = kang_ch_plot[frac < 0.035],
            aes(x = x_callout, y = y_callout, label = lbl, fill = channel_category),
            color = "#15151A", size = 2.0, fontface = "bold",
            label.size = 0.18, label.padding = unit(0.07, "lines"),
            min.segment.length = 0, segment.color = "#606876",
            segment.size = 0.18, seed = 42, show.legend = FALSE) +
  annotate("text", x = 0, y = 0.08, label = scales::comma(total_feat),
           size = 4.5, fontface = "bold", color = "#111827") +
  annotate("text", x = 0, y = -0.12, label = "active\ngenes",
           size = 2.25, lineheight = 0.78, color = "#4B5563") +
  scale_fill_manual(values = KANG_CH_COLS) +
  coord_fixed(xlim = c(-1.13, 1.13), ylim = c(-1.13, 1.13), clip = "off") +
  labs(title = "**A.**",
        x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(plot.title = ggtext::element_markdown(size = 11, face = "bold",
                                               hjust = 0, lineheight = 1.0),
        plot.margin = margin(4, 8, 2, 4))

# --- Panel B: DE-only GO BP enrichment ----------------------------------
kang_en_bp <- kang_en[channel_category == "DE-only" & database == "GO" &
                        ontology == "BP" & !is.na(p.adjust)]
kang_en_bp <- kang_en_bp[order(p.adjust)][seq_len(min(8L, .N))]
kang_en_bp[, desc_wrap := stringr::str_wrap(Description, width = 34)]
kang_en_bp[, desc_wrap := factor(desc_wrap, levels = rev(desc_wrap))]
kang_en_bp[, neglog10p := -log10(pmax(p.adjust, 1e-300))]
kang_en_bp[, q_lbl := sprintf("q=%s", formatC(p.adjust, format = "e", digits = 1))]

p8b <- ggplot(kang_en_bp, aes(x = neglog10p, y = desc_wrap)) +
  geom_segment(aes(x = 0, xend = neglog10p, yend = desc_wrap),
                color = "#6A8FBF", linewidth = 0.4) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
              color = "grey55", linewidth = 0.3) +
  geom_point(aes(size = Count),
              shape = 21, fill = "#355C8A", color = "grey20",
              stroke = 0.3, alpha = 0.92) +
  geom_text(aes(label = q_lbl), hjust = -0.12, size = 2.2, color = "grey25") +
  scale_size_continuous(range = c(2.4, 6.4), name = "Gene count") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(title = "**B.**",
        x = expression(-log[10](adj.~p)), y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.y = element_text(size = 6.8),
         legend.position = "bottom",
         legend.key.size = unit(0.35, "cm"),
         plot.margin = margin(4, 6, 2, 4))

# --- Panel C: cell-type marker enrichment -------------------------------
kang_mk2 <- copy(kang_mk)[channel_category %in% c("DE-only", "DV-only",
                                                     "DG-only", "overlap") &
                              !is.na(p.adjust) & n_overlap > 0]
kang_mk2 <- kang_mk2[order(channel_category, -fold_enrichment, p.adjust)]
kang_mk2 <- kang_mk2[, head(.SD, 5L), by = channel_category]
kang_mk2[, neglog10p := -log10(pmax(p.adjust, 1e-300))]
kang_mk2[, cell_type_wrap := stringr::str_wrap(cell_type, width = 20)]
kang_mk2 <- kang_mk2[order(fold_enrichment)]
kang_mk2[, cell_type_wrap := factor(cell_type_wrap, levels = cell_type_wrap)]
kang_mk2[, n_lbl := sprintf("n=%d", n_overlap)]

p8c <- ggplot(kang_mk2, aes(x = fold_enrichment, y = cell_type_wrap)) +
  geom_segment(aes(x = 0, xend = fold_enrichment, yend = cell_type_wrap),
                linewidth = 0.4, alpha = 0.6, color = "#3B4CC0") +
  geom_point(aes(size = neglog10p),
              shape = 21, fill = "#3B4CC0", color = "grey20",
              stroke = 0.3, alpha = 0.92) +
  geom_text(aes(label = n_lbl), hjust = -0.15, size = 2.2, color = "grey25") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55",
              linewidth = 0.3) +
  scale_size_continuous(range = c(2.4, 6.2),
                         name = expression(-log[10]~q)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(title = "**C.**",
        x = "Fold enrichment", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(legend.position = "bottom",
         legend.key.size = unit(0.35, "cm"),
         axis.text.y = element_text(size = 7),
         plot.margin = margin(2, 6, 4, 4))

# --- Panel D: representative violins (adaptive by non-empty channel) ----
kang_gene_cap <- c("DE-only" = 2L, "DV-only" = 1L, "DG-only" = 1L, "overlap" = 2L)
kang_dist2 <- kang_dist[channel_category %in% c("DE-only", "DV-only",
                                                   "DG-only", "overlap") &
                          is.finite(expression)]
kang_dist2 <- kang_dist2[rank_within_category <= kang_gene_cap[channel_category]]
kang_dist2[, channel_category := factor(channel_category,
  levels = c("DE-only", "DV-only", "DG-only", "overlap"))]
kang_dist2 <- kang_dist2[order(channel_category, rank_within_category)]
# Truncate long ENSG-embedded labels and rebuild compact facet strips.
kang_dist2[, gene_label := sub("_ENSG[0-9]+$", "", gene_label)]
kang_dist2[, gene_label := sub("^([^|]+ \\| )?([^_]+)_.*$", "\\1\\2", gene_label)]
kang_channel_short <- c("DE-only" = "DE", "DV-only" = "DV",
                        "DG-only" = "DG", "overlap" = "OL")
kang_dist2[, facet_lbl := paste(kang_channel_short[channel_category], gene_label, sep = " | ")]
kang_dist2[, facet_lbl := factor(facet_lbl, levels = unique(facet_lbl))]
grp_lvls <- unique(as.character(kang_dist2$group))
grp_cols <- c(ctrl = "#4C72B0", stim = "#C44E52")[grp_lvls]
if (any(is.na(grp_cols))) grp_cols[is.na(grp_cols)] <- "grey70"

p8d <- ggplot(kang_dist2, aes(x = group, y = expression, fill = group)) +
  ggdist::stat_halfeye(adjust = 0.8, width = 0.5, .width = 0,
                         justification = -0.15, point_colour = NA,
                         alpha = 0.5) +
  geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 0.9,
                color = "grey25", linewidth = 0.3) +
  facet_wrap(~facet_lbl, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(title = "**D.**",
        x = NULL, y = "Expression") +
  THEME_JBHI(base_size = 10) +
  theme(strip.text = element_text(size = 6.8, face = "bold", lineheight = 0.9),
         axis.text = element_text(size = 7),
         panel.spacing = unit(0.28, "lines"),
         plot.margin = margin(4, 4, 4, 4))

top_8 <- cowplot::plot_grid(p8a, p8b, ncol = 2,
                              rel_widths = c(0.9, 1.35), align = "h", axis = "tb")
bottom_8 <- cowplot::plot_grid(p8c, p8d, ncol = 2,
                                 rel_widths = c(0.82, 1.43), align = "h", axis = "tb")
fig8 <- cowplot::plot_grid(top_8, bottom_8,
                             ncol = 1, rel_heights = c(0.92, 1.04))
SAVE_FIG_V5("Fig8_kang18_merged", fig8, 7.5, 8.7)

# =========================================================================
# FigS1 — Manhattan-style channel decomposition (polished)
# =========================================================================
cat("  Fig S1 (v5): Manhattan...\n")
kirc_dt <- copy(kirc_dt_all)
kirc_dt[, ch_num := as.integer(channel)]
set.seed(12345)
kirc_dt[, x_jit := ch_num + runif(.N, -0.36, 0.36)]
y_cap <- quantile(abs(kirc_dt$log2FC), 0.995, na.rm = TRUE)
kirc_dt[, y_show := pmin(pmax(log2FC, -y_cap), y_cap)]

# Label selection: reduce overlap by selecting genes that are spatially separated
kirc_dt[, dist_score := abs(log2FC) * (-log10(pmax(padj_de, 1e-300)))]
nlbl_per_ch <- c("DE-only" = 6L, "DV-only" = 6L, "DG-only" = 4L,
                   "DE & DV" = 5L, "DD-only (not DE)" = 5L,
                   "Not significant" = 0L)
top_genes <- unlist(lapply(names(nlbl_per_ch), function(ch) {
  sub <- kirc_dt[channel == ch]
  if (nrow(sub) == 0L) return(character())
  head(sub[order(-dist_score)], nlbl_per_ch[[ch]])$gene
}))
top_genes <- unique(top_genes[!is.na(top_genes)])

CH_BG_V5 <- colorspace::lighten(CH_COLS_V5, 0.56)
bg_df <- data.table(
  xmin = seq_along(CH_LEVELS) - 0.5, xmax = seq_along(CH_LEVELS) + 0.5,
  channel = factor(CH_LEVELS, levels = CH_LEVELS))
ch_counts_s1 <- merge(
  data.table(channel = factor(CH_LEVELS, levels = CH_LEVELS)),
  kirc_dt[, .N, by = channel],
  by = "channel", all.x = TRUE)
ch_counts_s1[is.na(N), N := 0L]
ch_counts_s1[, label := paste0(channel, "\n(n = ",
                                  scales::comma(N), ")")]
lfc_lim <- max(abs(kirc_dt$y_show), na.rm = TRUE)

figS1 <- ggplot(kirc_dt, aes(x = x_jit, y = y_show)) +
  geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax,
                                ymin = -y_cap * 1.05,
                                ymax = y_cap * 1.05, fill = channel),
             inherit.aes = FALSE, show.legend = FALSE, alpha = 1) +
  geom_violin(data = kirc_dt[channel != "Not significant"],
              aes(x = ch_num, y = y_show, group = channel, fill = channel),
              inherit.aes = FALSE, width = 0.88, alpha = 0.42,
              color = "white", linewidth = 0.36, trim = TRUE) +
  geom_hline(yintercept = 0, color = "grey25", linewidth = 0.3) +
  rasterise(geom_point(aes(color = y_show),
                          size = 0.24, alpha = 0.34, shape = 16),
             dpi = 300) +
  geom_point(data = kirc_dt[gene %in% top_genes], aes(color = y_show),
              size = 1.8, alpha = 0.95, shape = 16) +
  geom_point(data = kirc_dt[gene %in% top_genes], color = "black",
              size = 2.2, alpha = 0.8, shape = 1, stroke = 0.55) +
  ggrepel::geom_text_repel(
     data = kirc_dt[gene %in% top_genes],
     aes(label = gene), size = 2.5, color = "black",
    bg.color = "white", bg.r = 0.16,
    max.overlaps = Inf, segment.size = 0.22, segment.color = "grey30",
    box.padding = 0.35, point.padding = 0.2, force = 20, force_pull = 0.15,
    seed = 42, min.segment.length = 0,
    segment.curvature = -0.2, segment.ncp = 3,
    nudge_y = fifelse(kirc_dt[gene %in% top_genes]$y_show >= 0,
                       y_cap * 0.20, -y_cap * 0.20)) +
  scale_color_gradientn(
    name = expression(log[2]*"FC"),
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-lfc_lim, lfc_lim), oob = scales::squish,
    breaks = c(-lfc_lim, 0, lfc_lim),
    labels = function(x) sprintf("%.1f", x),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
      barwidth = unit(0.5, "cm"), barheight = unit(3.2, "cm"),
      title.hjust = 0.5, title.position = "top")) +
  scale_fill_manual(values = CH_BG_V5) +
  scale_x_continuous(breaks = seq_along(CH_LEVELS),
                      labels = ch_counts_s1$label,
                      expand = expansion(mult = 0.01)) +
  coord_cartesian(ylim = c(-y_cap * 1.05, y_cap * 1.05), clip = "on") +
  labs(title = NULL,
        x = NULL,
        y = expression(Average ~ log[2] ~ "fold change")) +
  THEME_JBHI(base_size = 10) +
  theme(panel.grid = element_blank(),
         axis.text.x = element_text(face = "bold", size = 9, color = "black",
                                      margin = margin(t = 5), lineheight = 1.05),
         axis.ticks.x = element_blank(),
         panel.border = element_rect(color = "#15151A", fill = NA,
                                     linewidth = 0.35),
         legend.position = "right")
SAVE_FIG_V5("FigS1_tcga_kirc_manhattan", figS1, 14, 6.5)

# =========================================================================
# FigS2 — gene halfeye with |log2FC|/|log2CV2| annotations
# =========================================================================
cat("  Fig S2 (v5): gene boxplots...\n")
kirc_raw <- kirc_expr; counts_kirc2 <- counts_kirc; group_kirc2 <- group_kirc
de_cand2 <- kirc_dt_all[cat_de_only == TRUE & padj_de < 1e-6 &
                          padj_dv > 0.1][order(-abs(log2FC))]
de_top3 <- head(de_cand2, 3)$gene
dv_cand2 <- kirc_dt_all[cat_dv_only == TRUE & padj_de > 0.3 &
                          abs(log2FC) < 0.5][order(padj_dv)]
dv_top3 <- head(dv_cand2, 3)$gene
rep_genes2 <- c(de_top3, dv_top3)

expr_mat2 <- log2(nc_b[rep_genes2, , drop = FALSE])
gene_meta2 <- kirc_dt_all[gene %in% rep_genes2]
gene_meta2[, channel := fifelse(cat_de_only == TRUE, "DE-only", "DV-only")]
gene_meta2[, facet_lbl := paste0(gene, " (", channel, ")")]
gene_meta2[, lbl := sprintf("|log2FC|=%.2f \u00b7 |\u0394log2CV\u00b2|=%.2f",
                               abs(log2FC), abs(dv_log2_cv2_ratio))]
gene_meta2[, gene_order := match(gene, rep_genes2)]
expr_long2 <- data.table(
  gene = rep(rep_genes2, each = ncol(expr_mat2)),
  expr = as.vector(t(expr_mat2)),
  group = rep(as.character(group_kirc2), length(rep_genes2)))
expr_long2 <- merge(expr_long2,
                     gene_meta2[, .(gene, facet_lbl, lbl, gene_order)],
                     by = "gene")
expr_long2[, facet_lbl := reorder(facet_lbl, gene_order)]

figS2 <- ggplot(expr_long2, aes(x = group, y = expr, fill = group)) +
  ggdist::stat_halfeye(adjust = 0.8, width = 0.5, .width = 0,
                         justification = -0.18, point_colour = NA,
                         alpha = 0.55) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85,
                color = "grey25", linewidth = 0.3) +
  ggbeeswarm::geom_quasirandom(width = 0.08, size = 0.3, alpha = 0.25,
                                shape = 16) +
  geom_text(data = unique(gene_meta2[, .(facet_lbl, lbl)]),
             aes(x = 1.5, y = Inf, label = lbl),
             inherit.aes = FALSE, vjust = 1.4, size = 2.2,
             color = "grey25", fontface = "italic") +
  facet_wrap(~facet_lbl, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c(normal = "#2F66FF", tumor = "#FF375F"),
                      name = "Group") +
  labs(title = NULL,
        x = "Condition", y = expression(log[2]~expression)) +
  THEME_JBHI(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"),
         legend.position = "bottom")
SAVE_FIG_V5("FigS2_tcga_kirc_gene_boxplots", figS2, 13, 6.5)

# =========================================================================
# FigS3 — Pan-cancer channels (stacked bar + focus lollipop)
# =========================================================================
cat("  Fig S3 (v5): Pan-cancer...\n")
pan_ch <- fread(file = file.path(CONFIG$IN_DIR, "tcga_pan_cancer_sgcb_channels.csv"))
pan_l <- melt(pan_ch, id.vars = c("cancer", "n_genes"),
               measure.vars = c("de_only", "dv_only", "both_de_dv", "dd_only"),
               variable.name = "channel", value.name = "count")
PAN_CH <- c("DE-only", "DV-only", "DE & DV", "DD-only (not DE)")
pan_l[, channel := factor(channel,
  levels = c("de_only", "dv_only", "both_de_dv", "dd_only"),
  labels = PAN_CH)]
c_ord <- pan_ch[order(-dv_only), cancer]
pan_l[, cancer := factor(cancer, levels = c_ord)]
pan_ch[, cancer := factor(cancer, levels = c_ord)]
pan_l[, proportion := count / sum(count), by = cancer]
pan_l[, plabel := scales::percent(proportion, accuracy = 0.1)]
pan_l[, txt_col := fifelse(proportion >= 0.38, "white", "#15151A")]

pS3a <- ggplot(pan_l, aes(x = cancer, y = channel)) +
  geom_tile(fill = "#F3F6FA", color = "white", linewidth = 0.45) +
  geom_point(aes(size = proportion, fill = channel),
             shape = 21, color = "white", stroke = 0.5, alpha = 0.96) +
  geom_text(aes(label = plabel, color = txt_col),
            size = 2.35, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_manual(values = CH_COLS_V5[PAN_CH], name = "Channel") +
  scale_size_area(max_size = 15,
                  breaks = c(0.10, 0.25, 0.50),
                  labels = scales::percent_format(accuracy = 1),
                  name = "Proportion") +
  labs(title = "**A.**",
        x = NULL, y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "#15151A", fill = NA,
                                    linewidth = 0.35))

pan_foc <- pan_l[channel %in% c("DV-only", "DD-only (not DE)")]
pS3b <- ggplot(pan_foc, aes(x = cancer, y = count, color = channel)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  geom_segment(aes(xend = cancer, y = 0, yend = count),
                linewidth = 0.35, alpha = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = scales::comma(count)),
             vjust = -0.8, size = 2.2, show.legend = FALSE) +
  ggh4x::facet_grid2(rows = vars(channel), scales = "free_y",
                       independent = "y", switch = "y") +
  scale_color_manual(values = CH_COLS_V5[c("DV-only", "DD-only (not DE)")],
                      guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "**B.**",
        x = NULL, y = "Genes") +
  THEME_JBHI(base_size = 10) +
  theme(strip.placement = "outside",
         strip.background = element_rect(fill = "#EEF2F8", color = NA),
         strip.text.y.left = element_text(angle = 0, face = "bold", size = 7),
         axis.text.x = element_text(angle = 30, hjust = 1))

figS3 <- (pS3a / pS3b) + plot_layout(heights = c(1.2, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"))
SAVE_FIG_V5("FigS3_pan_cancer_channels", figS3, 10, 9)

# =========================================================================
# FigS4 — null calibration and robustness checks (2x2 compact)
# =========================================================================
cat("  Fig S4 (v5): robustness checks (2x2)...\n")
REV_DIR <- CONFIG$REVIEW_DIR

# A: Null FPR by channel
rv_null <- fread(file = file.path(REV_DIR, "A_null_channel_fpr_summary.csv"))
rv_null[, channel := factor(channel,
  levels = c("DE", "DV", "DG_gamma", "DD", "SGCB_Score"))]
rv_null[, channel_disp := factor(as.character(channel),
  levels = c("DE", "DV", "DG_gamma", "DD", "SGCB_Score"),
  labels = c("DE", "DV", "DG-gamma", "DD", "Score"))]
rv_null[, lbl := sprintf("%.4f", mean_fpr)]
pS4a <- ggplot(rv_null, aes(x = mean_fpr, y = channel_disp)) +
  annotate("rect", xmin = 0, xmax = 0.05, ymin = -Inf, ymax = Inf,
           fill = "#EEF5FF", alpha = 0.65) +
  annotate("rect", xmin = 0.05, xmax = 0.06, ymin = -Inf, ymax = Inf,
           fill = "#FFF1F1", alpha = 0.75) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#B22222",
              linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = mean_fpr, y = channel_disp, yend = channel_disp),
                color = "#6F7D95", linewidth = 1.0, lineend = "round") +
  geom_errorbar(aes(xmin = pmax(mean_fpr - sd_fpr, 0), xmax = q95_fpr),
                 orientation = "y", width = 0.15,
                 linewidth = 0.35, color = "#5B677A") +
  geom_point(aes(fill = channel_disp), size = 3.4, shape = 21,
             color = "white", stroke = 0.7) +
  geom_text(aes(label = lbl), hjust = -0.35, size = 2.25, color = "grey25") +
  annotate("text", x = 0.049, y = 5.45, label = "5% target",
            hjust = 1, size = 2.2, color = "#B22222", fontface = "italic") +
  scale_fill_manual(values = c("DE" = "#3157FF", "DV" = "#25C7AE",
                                "DG-gamma" = "#7A4CFF", "DD" = "#FF8C3A",
                                "Score" = "#0F172A"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 0.06), breaks = c(0, 0.01, 0.025, 0.05),
                      labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "**A.**",
        x = "Empirical FPR", y = NULL) +
  THEME_JBHI(base_size = 10)

# B: Permutation baseline (observed vs permuted)
rv_perm <- fread(file = file.path(REV_DIR, "C_tcga_permutation_test_summary.csv"))
rv_perm[, endpoint := factor(endpoint,
  levels = c("DV_only", "DD_only", "shape_proxy"),
  labels = c("DV-only", "DD-only", "Shape proxy"))]
rv_perm[, observed_n := as.numeric(observed_n)]
rv_perm[, perm_plot := perm_mean + 1]
rv_perm[, obs_plot := observed_n + 1]
rv_perm[, obs_lab := scales::comma(observed_n, accuracy = 1)]
rv_perm[, perm_lab := sprintf("%.1f", perm_mean)]
rv_perm[, ratio := observed_n / pmax(perm_mean, .Machine$double.eps)]
rv_perm[, ratio_lab := fifelse(ratio < 10,
  sprintf("%.1fx", ratio),
  sprintf("%sx", scales::comma(ratio, accuracy = 1)))]
pS4b <- ggplot(rv_perm, aes(y = endpoint)) +
  geom_segment(aes(x = perm_plot, xend = obs_plot,
                   yend = endpoint),
               linewidth = 1.35, color = "#91A7FF",
               alpha = 0.72, lineend = "round") +
  geom_point(aes(x = perm_plot), size = 3.2, shape = 21,
             fill = "white", color = "#70819B", stroke = 0.9) +
  geom_point(aes(x = obs_plot), size = 4.0, shape = 21,
             fill = "#3157FF", color = "white", stroke = 0.8) +
  geom_text(aes(x = perm_plot, label = perm_lab),
            hjust = 1.35, size = 2.15, color = "#5B677A") +
  geom_text(aes(x = obs_plot, label = obs_lab),
            hjust = -0.25, size = 2.2, color = "#1F2A44",
            fontface = "bold") +
  geom_text(aes(x = sqrt(perm_plot * obs_plot), label = ratio_lab),
            vjust = -0.9, size = 2.05, color = "#3157FF",
            fontface = "bold") +
  scale_x_log10(labels = scales::comma,
                breaks = c(1, 10, 100, 1000, 3000),
                expand = expansion(mult = c(0.08, 0.2))) +
  labs(title = "**B.**",
        x = "Genes + 1 (log scale)", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

# C: Dose-response monotonicity
rv_mono <- fread(file = file.path(REV_DIR, "B_dose_response_monotonicity.csv"))
rv_mono[, endpoint_label := fifelse(endpoint == "DV_channel", "DV ch.",
  fifelse(endpoint == "DV_only_not_DE", "DV-only",
  fifelse(endpoint == "DD_channel" & target == "truth_DV", "DD vs DV",
  fifelse(endpoint == "DD_channel" & target == "truth_DV_or_DD", "DD vs DV/DD",
  fifelse(endpoint == "DD_only_not_DE", "DD-only", "DD shape")))))]
rv_mono[, endpoint_label := factor(endpoint_label,
  levels = c("DV ch.", "DV-only", "DD vs DV", "DD vs DV/DD",
              "DD-only", "DD shape"))]
rv_ml <- melt(rv_mono, id.vars = c("dose_type", "endpoint", "target",
                                     "endpoint_label"),
               measure.vars = c("spearman_rho_recall", "pearson_r_recall"),
               variable.name = "metric", value.name = "value")
rv_ml[, metric := factor(metric, levels = c("spearman_rho_recall",
                                               "pearson_r_recall"),
  labels = c("Spearman", "Pearson"))]
rv_ml[, value_lab := sprintf("%.2f", value)]

pS4c <- ggplot(rv_ml, aes(x = value, y = endpoint_label,
                           color = dose_type, shape = metric)) +
  annotate("rect", xmin = 0.8, xmax = 1.08, ymin = -Inf, ymax = Inf,
           fill = "#EAF8F4", alpha = 0.8) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_vline(xintercept = 0.8, linetype = "dashed",
              color = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = 0.8, xend = value, yend = endpoint_label),
               position = position_dodge(width = 0.45),
               linewidth = 0.75, alpha = 0.42, lineend = "round") +
  geom_point(position = position_dodge(width = 0.45), size = 2.5,
              stroke = 0.8) +
  geom_text(aes(label = value_lab),
             position = position_dodge(width = 0.45),
             hjust = -0.35, size = 2.05, color = "grey20") +
  scale_color_manual(values = c(DV = "#2F5597", DD = "#2A9D8F"),
                      name = "Injected channel") +
  scale_shape_manual(values = c(Spearman = 16, Pearson = 1), name = NULL) +
  scale_x_continuous(limits = c(-0.05, 1.08),
                      breaks = c(0, 0.5, 0.8, 1.0)) +
  labs(title = "**C.**",
        x = "Recall correlation across dose levels", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.y = element_text(size = 7),
         legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

# D: Subtype imbalance (balanced vs imbalanced)
rv_conf <- fread(file = file.path(REV_DIR, "D_shape_variance_confound_contrast.csv"))
rv_conf <- rv_conf[metric %in% c("mean_dv_fpr", "mean_dd_fpr",
                                     "mean_n_dv_only", "mean_n_dd_only",
                                     "mean_dv_confounded_capture",
                                     "mean_dd_confounded_capture")]
rv_conf[, metric_label := factor(metric,
  levels = c("mean_dv_fpr", "mean_dd_fpr",
              "mean_dv_confounded_capture", "mean_dd_confounded_capture",
              "mean_n_dv_only", "mean_n_dd_only"),
  labels = c("DV FPR", "DD FPR",
              "DV confound", "DD confound",
              "DV-only genes", "DD-only genes"))]
rv_conf[, metric_group := fifelse(metric %in% c("mean_n_dv_only",
                                                   "mean_n_dd_only"),
                                     "Gene counts", "Proportions")]
rv_long <- melt(rv_conf, id.vars = c("metric", "metric_label", "metric_group"),
                 measure.vars = c("balanced", "imbalanced"),
                 variable.name = "condition", value.name = "value")
rv_long[, condition := factor(condition, levels = c("balanced", "imbalanced"),
  labels = c("Balanced", "Imbalanced"))]
rv_long[, value_lab := fifelse(metric_group == "Gene counts",
  scales::comma(value, accuracy = 0.1),
  scales::percent(value, accuracy = 0.01))]
rv_long[, value_lab := fifelse(value == 0, "0", value_lab)]

pS4d <- ggplot(rv_long, aes(x = value, y = metric_label,
                             color = condition, shape = condition)) +
  geom_segment(aes(x = 0, xend = value, yend = metric_label),
               position = position_dodge(width = 0.55),
               linewidth = 0.8, alpha = 0.42, lineend = "round") +
  geom_point(position = position_dodge(width = 0.55), size = 2.5,
              stroke = 0.8) +
  geom_text(aes(label = value_lab),
             position = position_dodge(width = 0.55),
             hjust = -0.25, size = 2.05, color = "grey15") +
  ggh4x::facet_grid2(cols = vars(metric_group), scales = "free_x",
                       independent = "x", switch = "x") +
  scale_color_manual(values = c(Balanced = "#4C72B0", Imbalanced = "#D4A017"),
                      name = "Subtype mix") +
  scale_shape_manual(values = c(Balanced = 16, Imbalanced = 17),
                      name = "Subtype mix") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "**D.**",
        x = NULL, y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(strip.placement = "outside",
         strip.background = element_rect(fill = "grey96", color = NA),
         axis.text.y = element_text(size = 7),
         legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

figS4 <- (pS4a | pS4b) / (pS4c | pS4d) +
  plot_layout(widths = c(1.03, 1), heights = c(1, 1.08))
SAVE_FIG_V5("FigS4_reviewer_checks", figS4, 13, 10)

cat("  37c_v5: Biological + supplementary figures done.\n")
