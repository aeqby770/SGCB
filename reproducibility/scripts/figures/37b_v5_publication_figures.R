# =========================================================================
# 37b_v5 — SGCB JBHI paper Fig 2–6 + Fig7 SEQC (benchmark figures)
# 数据由 37_v5_run_all.R 入口加载
# 关键升级:
#   - Fig3 F1 用 ComplexHeatmap + 列 split by effect size, 行注释 mean F1
#   - Fig4 practical 合并成单张 ComplexHeatmap (classic datasets heatmap
#     + F1 across n 曲线 + FDR annotation), 解决 lollipop imbalance 问题
#   - Fig6 multiomics 改为 3 panel: pseudobulk-heatmap / prot-heatmap / spike-in
#     每个 heatmap 右侧 column annotation = null FDR barplot
#   - 统一 ggsci NPG 配色 + ggpubr theme_pubr 骨架
# =========================================================================

source(file.path(PROJECT_ROOT, "scripts", "figures", "37_theme_v5.R"),
       local = FALSE, encoding = "UTF-8")

# =========================================================================
# Fig 2 — Null Calibration (hist + QQ + FDR band)
# =========================================================================
cat("  Fig 2 (v5): Null Calibration...\n")
null_n10 <- copy(combined_null_pvalues)
null_n10[, method := factor(method, levels = METHOD_ORDER_V3)]

# Panel A — p-value density fingerprint. A compact heat-strip makes the
# anti-conservative edge pile-ups visually immediate without repeating axes.
null_bins <- null_n10[, .N, by = .(
  method,
  bin = pmin(49L, as.integer(floor(pvalue * 50)))
)]
null_bins[, total := sum(N), by = method]
null_bins[, density := N / (total * 0.02)]
null_bins[, bin_mid := (bin + 0.5) / 50]
null_bins[, method := factor(method, levels = rev(METHOD_ORDER_V3))]

ks_dt <- null_n10[, .(ks = round(ks.test(pvalue, "punif")$statistic, 3)),
                    by = method]
ks_dt[, method := factor(method, levels = rev(METHOD_ORDER_V3))]
ks_dt[, label := sprintf("KS %.3f", ks)]

p2a <- ggplot(null_bins, aes(x = bin_mid, y = method, fill = density)) +
  geom_tile(width = 0.020, height = 0.78, color = "white", linewidth = 0.08) +
  geom_vline(xintercept = c(0.05, 0.95), color = "white", linewidth = 0.28,
             linetype = "dashed") +
  geom_text(data = ks_dt, aes(x = 1.04, y = method, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.25, color = "#1B1B1F",
            fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#F8FBFF", "#C7F2EA", "#34D6C8", "#3157E8", "#281447"),
    limits = c(0, quantile(null_bins$density, 0.985, na.rm = TRUE)),
    oob = scales::squish,
    name = "Density") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.18),
                     expand = expansion(0, 0)) +
  labs(title = "**A.** Null p-value fingerprint",
        x = "p-value", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8.4),
        legend.position = "bottom",
        legend.key.width = unit(1.05, "cm"),
        plot.margin = margin(4, 6, 4, 8))

# Panel B — Q-Q deviation fingerprint. Blue means conservative upward p-values,
# red means left-tail inflation; SGCB should sit close to the zero contour.
qq_dev <- null_n10[, {
  n <- .N
  theoretical <- ppoints(n)
  observed <- sort(pvalue)
  data.table(
    qbin = pmin(99L, as.integer(floor(theoretical * 100))),
    deviation = observed - theoretical
  )[, .(dev = mean(deviation)), by = qbin]
}, by = method]
qq_dev[, q_mid := (qbin + 0.5) / 100]
qq_dev[, method := factor(method, levels = rev(METHOD_ORDER_V3))]

p2b <- ggplot(qq_dev, aes(x = q_mid, y = method, fill = dev)) +
  geom_tile(width = 0.010, height = 0.78, color = NA) +
  geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "white",
             linewidth = 0.18, alpha = 0.75) +
  scale_fill_gradient2(
    low = "#E9435D", mid = "#F9FAFC", high = "#2156D9",
    midpoint = 0, limits = c(-0.35, 0.35), oob = scales::squish,
    name = "Observed - expected") +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = expansion(0, 0)) +
  labs(title = "**B.** Q-Q deviation field",
        x = "Theoretical quantile", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8.4),
        legend.position = "bottom",
        legend.key.width = unit(1.05, "cm"),
        plot.margin = margin(4, 6, 4, 8))

# Panel C — FDR calibration band
fdr_cal8 <- melt(combined_sim[, .(method, actual_FDR_fdr1, actual_FDR_fdr5,
                                    actual_FDR_fdr10)],
                  id.vars = "method",
                  measure.vars = c("actual_FDR_fdr1", "actual_FDR_fdr5",
                                    "actual_FDR_fdr10"),
                  variable.name = "threshold", value.name = "observed_fdr")
fdr_cal8[, nominal := fifelse(threshold == "actual_FDR_fdr1", 0.01,
                       fifelse(threshold == "actual_FDR_fdr5", 0.05, 0.10))]
fdr_summary <- fdr_cal8[, .(mean_fdr = mean(observed_fdr, na.rm = TRUE),
                             sd_fdr = sd(observed_fdr, na.rm = TRUE)),
                          by = .(method, nominal)]
fdr_summary[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
fdr_summary[is.na(sd_fdr), sd_fdr := 0]
calibrated <- c("SGCB", "limma", "glmGamPoi")
fdr_summary[, group := fifelse(method %in% calibrated, "Calibrated",
                                "Anti-conservative")]

p2c <- ggplot(fdr_summary, aes(x = nominal, y = mean_fdr, color = method)) +
  annotate("rect", xmin = 0, xmax = 0.105, ymin = 0, ymax = 0.05,
           fill = "#B9F6E3", alpha = 0.18) +
  annotate("rect", xmin = 0, xmax = 0.105, ymin = 0.05, ymax = Inf,
           fill = "#FF6D6D", alpha = 0.08) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35",
               linewidth = 0.35) +
  geom_ribbon(aes(ymin = pmax(mean_fdr - sd_fdr, 0),
                   ymax = mean_fdr + sd_fdr, fill = method, group = method),
               alpha = 0.10, color = NA) +
  geom_line(aes(linewidth = method == "SGCB"), alpha = 0.93) +
  geom_point(aes(size = method == "SGCB"), shape = 21, fill = "white",
             stroke = 0.9) +
  facet_wrap(~group, scales = "free_y") +
  scale_color_manual(values = METHOD_COLORS_V3, name = "Method") +
  scale_fill_manual(values = METHOD_COLORS_V3, guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.42, `TRUE` = 0.95), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.0, `TRUE` = 3.1), guide = "none") +
  scale_x_continuous(breaks = c(0.01, 0.05, 0.10),
                      labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "**C.** Calibration operating zone",
        x = "Nominal FDR", y = "Observed FDR") +
  THEME_JBHI(base_size = 10) +
  theme(legend.key.size = unit(0.35, "cm"),
         strip.text = element_text(face = "bold"),
         plot.margin = margin(4, 8, 4, 8))

fig2 <- cowplot::plot_grid(
  p2a, p2b, p2c,
  ncol = 1, rel_heights = c(0.86, 0.86, 1.0))
SAVE_FIG_V5("Fig2_calibration", fig2, 7.5, 8.2)

# =========================================================================
# Fig 3 — Simulation Power (ComplexHeatmap F1 + FDR + trajectory)
# =========================================================================
cat("  Fig 3 (v5): Simulation Power...\n")
sim_fig <- copy(combined_sim)
sim_fig[, n_per_group := as.integer(sub("sim_n([0-9]+)_.*", "\\1", dataset))]
sim_fig[, effect_size := as.numeric(sub(".*_ef([0-9.]+)$", "\\1", dataset))]
sim_agg <- sim_fig[, .(mean_F1  = mean(F1_fdr5, na.rm = TRUE),
                         mean_FDR = mean(actual_FDR_fdr5, na.rm = TRUE),
                         sd_FDR   = sd(actual_FDR_fdr5, na.rm = TRUE)),
                     by = .(method, n_per_group, effect_size)]
sim_agg[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
sim_agg[is.na(sd_FDR), sd_FDR := 0]

# Panel A — F1 heatmap (ggplot2 geom_tile + side-bar, NPG-styled)
sim_long <- copy(sim_agg)
sim_long[, ef_lab := paste0("Effect = ", effect_size, "\u00d7")]
# method order by mean F1 at (n=10, effect=2x)
ord_tbl <- sim_long[n_per_group == 10 & effect_size == 2,
                      .(method, mean_F1)][order(-mean_F1)]
meth_ord <- as.character(ord_tbl$method)
sim_long[, method := factor(method, levels = rev(meth_ord))]
sim_long[, F1_txt_col := ifelse(mean_F1 > 0.45, "white", "grey15")]
sgcb_sim_tiles <- sim_long[as.character(method) == "SGCB"]

p3a_heat <- ggplot(sim_long, aes(x = factor(n_per_group), y = method,
                                    fill = mean_F1)) +
  geom_tile(color = "white", linewidth = 0.38) +
  geom_tile(data = sgcb_sim_tiles, aes(x = factor(n_per_group), y = method),
            inherit.aes = FALSE, fill = NA, color = "#FF3B30",
            linewidth = 0.68) +
  geom_text(aes(label = sprintf("%.2f", mean_F1), color = F1_txt_col),
              size = 2.6, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(colors = c("#FAFCFF", "#D9FFF7", "#52E4CF", "#2A70FF", "#211442"),
                        limits = c(0, 1), name = "F1",
                        breaks = c(0, 0.5, 1),
                        guide = guide_colorbar(barwidth = 10, barheight = 0.4,
                          title.position = "left", title.vjust = 1)) +
  facet_wrap(~ef_lab, ncol = 2) +
  scale_x_discrete(expand = expansion(0, 0)) +
  scale_y_discrete(expand = expansion(0, 0)) +
  labs(title = "**A.**",
        x = "Samples per group", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(strip.text = element_text(size = 8.5, face = "bold"),
         legend.position = "bottom", legend.justification = "left",
         panel.border = element_rect(color = "grey30", fill = NA,
                                       linewidth = 0.3),
         axis.text.y = element_text(face = "bold", size = 8.5))

# Right-side mean F1 bar
row_ann_vals <- sim_agg[, .(mean_F1 = mean(mean_F1, na.rm = TRUE)),
                          by = method]
row_ann_vals[, method := factor(as.character(method),
                                  levels = rev(meth_ord))]
p3a_bar <- ggplot(row_ann_vals, aes(x = mean_F1, y = method, fill = method)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.2,
            show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", mean_F1)),
             hjust = -0.12, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = METHOD_COLORS_V3) +
  scale_x_continuous(limits = c(0, max(row_ann_vals$mean_F1) * 1.32),
                      expand = expansion(0, 0)) +
  labs(title = NULL, x = NULL, y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         plot.title = ggtext::element_markdown(size = 8.5, face = "bold",
                                                  hjust = 0.5, lineheight = 1.0),
         panel.border = element_rect(color = "grey30", fill = NA,
                                       linewidth = 0.3),
         plot.margin = margin(4, 8, 4, 2))
p3a <- cowplot::plot_grid(p3a_heat, p3a_bar, ncol = 2,
                            rel_widths = c(4.5, 1), align = "h",
                            axis = "tb")

# Panel B — FDR cap-at-25% with ggforce split axis
p3b_data <- copy(sim_agg)
p3b_data[, ef_lab := paste0("Effect = ", effect_size, "\u00d7")]
p3b_data[, fdr_disp := pmin(mean_FDR, 0.25)]
p3b_data[, capped := mean_FDR > 0.25]
p3b <- ggplot(p3b_data, aes(x = factor(n_per_group), y = fdr_disp,
                               color = method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
              color = FDR_LINE_COLOR, linewidth = 0.3) +
  geom_pointrange(aes(ymin = pmax(fdr_disp - sd_FDR, 0),
                       ymax = pmin(fdr_disp + sd_FDR, 0.25)),
                   position = position_dodge(width = 0.6),
                   size = 0.3, linewidth = 0.32) +
  geom_point(data = p3b_data[capped == TRUE],
              aes(x = factor(n_per_group), y = 0.245, color = method),
              shape = 17, size = 1.8, inherit.aes = FALSE) +
  facet_wrap(~ef_lab, ncol = 2) +
  scale_color_manual(values = METHOD_COLORS_V3, name = "Method") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                      breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25),
                      limits = c(0, 0.26), expand = expansion(0, 0)) +
  labs(title = "**B.**",
        x = "Samples per group", y = "Observed FDR") +
  THEME_JBHI(base_size = 10) + theme(legend.position = "right")

fig3 <- cowplot::plot_grid(
  p3a, p3b,
  ncol = 1, rel_heights = c(1.35, 1))
SAVE_FIG_V5("Fig3_simulation_power", fig3, 7.5, 8.3)

# =========================================================================
# Fig 4 — Practical operating range on classic datasets
# =========================================================================
cat("  Fig 4 (v5): Practical...\n")
classic8 <- copy(combined_classic)
classic8[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
classic_summary <- classic8[, .(n_sig = sum(n_sig_fdr5, na.rm = TRUE)),
                               by = .(method, dataset)]
classic_summary[, log_n := log10(pmax(n_sig, 0.1))]
classic_summary[n_sig == 0, log_n := 0]
classic_summary[, dataset_clean := clean_dataset_label(dataset)]
# Row order by mean yield
classic_mean <- classic_summary[, .(mean_yield = mean(n_sig, na.rm = TRUE)),
                                   by = method][order(-mean_yield)]
meth_ord_c <- as.character(classic_mean$method)
classic_summary[, method := factor(method, levels = rev(meth_ord_c))]
classic_summary[, txt_col := ifelse(log_n > 2.5, "white", "grey15")]
sgcb_classic_tiles <- classic_summary[as.character(method) == "SGCB"]

p4a_heat <- ggplot(classic_summary, aes(x = dataset_clean, y = method,
                                            fill = log_n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(data = sgcb_classic_tiles, aes(x = dataset_clean, y = method),
            inherit.aes = FALSE, fill = NA, color = "#FF3B30",
            linewidth = 0.68) +
  geom_text(aes(label = ifelse(n_sig == 0, "0", scales::comma(n_sig)),
                 color = txt_col),
              size = 2.4, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(colors = c("#FAFCFF", "#D9FFF7", "#52E4CF", "#2A70FF", "#211442"),
                        limits = c(0, max(classic_summary$log_n) * 1.01),
                        name = expression(log[10]*"(DE genes)"),
                        guide = guide_colorbar(barwidth = 0.4, barheight = 6,
                                                title.position = "top")) +
  scale_x_discrete(expand = expansion(0, 0)) +
  scale_y_discrete(expand = expansion(0, 0)) +
  labs(title = "**A.**",
        x = NULL, y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
         axis.text.y = element_text(face = "bold", size = 8.5),
         legend.position = "right",
         panel.border = element_rect(color = "grey30", fill = NA,
                                       linewidth = 0.3))

classic_mean[, method := factor(as.character(method),
                                   levels = rev(meth_ord_c))]
p4a_bar <- ggplot(classic_mean, aes(x = mean_yield, y = method,
                                       fill = method)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.2,
            show.legend = FALSE) +
  geom_text(aes(label = scales::comma(round(mean_yield))),
             hjust = -0.12, size = 2.3, fontface = "bold") +
  scale_fill_manual(values = METHOD_COLORS_V3) +
  scale_x_continuous(limits = c(0, max(classic_mean$mean_yield) * 1.4),
                      labels = scales::comma,
                      expand = expansion(0, 0)) +
  labs(title = NULL, x = NULL, y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         plot.title = ggtext::element_markdown(size = 8, face = "bold",
                                                  hjust = 0.5, lineheight = 1.0),
         panel.border = element_rect(color = "grey30", fill = NA,
                                       linewidth = 0.3),
         plot.margin = margin(4, 8, 4, 2))
p4a <- cowplot::plot_grid(p4a_heat, p4a_bar, ncol = 2,
                            rel_widths = c(5, 1), align = "h",
                            axis = "tb")

# Panel B — F1 across sample size (line, SGCB in bold red)
small_s <- copy(combined_sim)
small_s[, n_per_group := as.integer(sub("sim_n([0-9]+)_.*", "\\1", dataset))]
small_s <- small_s[, .(mean_F1 = mean(F1_fdr5, na.rm = TRUE),
                         sd_F1  = sd(F1_fdr5, na.rm = TRUE),
                         mean_FDR = mean(actual_FDR_fdr5, na.rm = TRUE)),
                     by = .(method, n_per_group)]
small_s[is.na(sd_F1), sd_F1 := 0]
small_s[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
small_s[, is_sgcb := method == "SGCB"]

p4b <- ggplot(small_s, aes(x = factor(n_per_group), y = mean_F1,
                              group = method, color = method)) +
  geom_ribbon(aes(ymin = pmax(mean_F1 - sd_F1, 0),
                   ymax = pmin(mean_F1 + sd_F1, 1),
                   fill = method), alpha = 0.10, color = NA) +
  geom_line(aes(linewidth = is_sgcb), alpha = 0.9) +
  geom_point(aes(size = is_sgcb), shape = 21, fill = "white", stroke = 0.9) +
  scale_color_manual(values = METHOD_COLORS_V3, name = "Method") +
  scale_fill_manual(values = METHOD_COLORS_V3, guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 1.1), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 3.2), guide = "none") +
  labs(title = "**B.**",
        x = "Samples per group", y = "Mean F1 \u00b1 SD") +
  THEME_JBHI(base_size = 10) + theme(legend.position = "right")

# Panel C — FDR heatmap across sample size
small_s_c <- copy(small_s)
small_s_c[, method := factor(as.character(method), levels = rev(meth_ord_c))]
small_s_c[, txt_col := ifelse(mean_FDR > 0.08, "white", "grey15")]

p4c <- ggplot(small_s_c, aes(x = factor(n_per_group), y = method,
                                fill = mean_FDR)) +
  geom_tile(color = "white", linewidth = 0.38) +
  geom_text(aes(label = sprintf("%.3f", mean_FDR), color = txt_col),
              size = 2.4, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(colors = c("#FAFCFF", "#B9F6E3", "#FFCA7A", "#FF5A5F", "#5B1022"),
                        values = scales::rescale(c(0, 0.05, 0.2, 0.7)),
                        limits = c(0, max(small_s_c$mean_FDR, na.rm = TRUE) * 1.01),
                        name = "Obs. FDR",
                        guide = guide_colorbar(barwidth = 0.4, barheight = 5,
                                                title.position = "top")) +
  labs(title = "**C.**",
        x = "Samples per group", y = NULL) +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
         legend.position = "right",
         panel.border = element_rect(color = "grey30", fill = NA,
                                       linewidth = 0.3))

p4bc <- cowplot::plot_grid(p4b, p4c, ncol = 2, rel_widths = c(1.4, 1),
                             align = "h", axis = "tb")
fig4 <- cowplot::plot_grid(
  p4a, p4bc,
  ncol = 1, rel_heights = c(1.3, 1))
SAVE_FIG_V5("Fig4_practical", fig4, 7.5, 9.3)

# =========================================================================
# Fig 6 — Multiomics: pseudobulk + proteomics + spike-in F1
# =========================================================================
cat("  Fig 6 (v5): Multiomics...\n")

# ---- helper: build (heatmap + null-FDR side-bar) for one modality -------
build_modality_panel <- function(de_dt, null_dt, strip_prefix_re,
                                    method_palette,
                                    panel_letter, modality_label) {
  de_dt  <- copy(de_dt)[!grepl("^dv_", method) & !is.na(method) & method != ""]
  null_dt<- copy(null_dt)[!grepl("^dv_", method) & !is.na(method) & method != ""]
  de_dt[, dataset_clean := clean_dataset_label(gsub(strip_prefix_re, "", dataset))]
  de_dt <- de_dt[, .(n_sig = sum(n_sig_fdr5, na.rm = TRUE)),
                   by = .(method, dataset_clean)]
  de_dt[, log_n := log10(pmax(n_sig, 0.1))]
  de_dt[n_sig == 0, log_n := 0]
  # row order by mean yield
  row_mean <- de_dt[, .(mean_yield = mean(n_sig, na.rm = TRUE)),
                      by = method][order(-mean_yield)]
  m_ord <- as.character(row_mean$method)
  de_dt[, method := factor(method, levels = rev(m_ord))]
  de_dt[, txt_col := ifelse(log_n > 2.5, "white", "grey15")]

  sgcb_tiles <- de_dt[as.character(method) == "SGCB"]
  p_heat <- ggplot(de_dt, aes(x = dataset_clean, y = method, fill = log_n)) +
    geom_tile(color = "white", linewidth = 0.42) +
    geom_tile(data = sgcb_tiles, aes(x = dataset_clean, y = method),
              inherit.aes = FALSE, fill = NA, color = "#FF3B30",
              linewidth = 0.72) +
    geom_text(aes(label = ifelse(n_sig == 0, "0", scales::comma(n_sig)),
                   color = txt_col),
               size = 2.25, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = c("#FAFCFF", "#D9FFF7", "#52E4CF", "#2A70FF", "#211442"),
      limits = c(0, max(de_dt$log_n, na.rm = TRUE) * 1.01),
      name = expression(log[10]*"(sig)"),
      guide = guide_colorbar(barwidth = 2.1, barheight = 0.22,
                               title.position = "top")) +
    scale_x_discrete(expand = expansion(0, 0)) +
    scale_y_discrete(expand = expansion(0, 0)) +
    labs(title = sprintf("**%s.**", panel_letter),
          x = NULL, y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
           axis.text.y = element_text(face = "bold", size = 8),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.title = element_text(size = 7.2, face = "bold"),
           legend.text = element_text(size = 6.8),
           panel.border = element_rect(color = "grey30", fill = NA,
                                         linewidth = 0.3),
           plot.margin = margin(4, 3, 2, 4))
  # Side bar: mean null FDR per method
  null_mean <- null_dt[, .(fdr = mean(mean_FDR05, na.rm = TRUE)),
                          by = method]
  null_mean <- null_mean[match(m_ord, method)]
  null_mean[, method := factor(m_ord, levels = rev(m_ord))]
  null_mean[, exceed := fdr > 0.05]
  p_bar <- ggplot(null_mean, aes(x = fdr, y = method,
                                    fill = exceed)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.25,
              show.legend = FALSE) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "#FF1744",
                linewidth = 0.5) +
    geom_point(aes(x = fdr), shape = 21, fill = "white", color = "#15151A",
               stroke = 0.32, size = 1.55, show.legend = FALSE) +
    scale_fill_manual(values = c(`FALSE` = "#2F66FF", `TRUE` = "#FF375F")) +
    scale_x_continuous(limits = c(0, max(c(null_mean$fdr, 0.06),
                                            na.rm = TRUE) * 1.55),
                        breaks = c(0, 0.05),
                        labels = scales::percent_format(accuracy = 0.1),
                        expand = expansion(0, 0)) +
    labs(title = NULL,
          x = NULL, y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
           plot.title = ggtext::element_markdown(size = 7.5, face = "bold",
                                                    hjust = 0.5, lineheight = 1.0),
           panel.border = element_rect(color = "grey30", fill = NA,
                                         linewidth = 0.3),
           plot.margin = margin(4, 8, 2, 2))
  cowplot::plot_grid(p_heat, p_bar, ncol = 2, rel_widths = c(5.0, 1.18),
                       align = "h", axis = "tb")
}

p6a <- build_modality_panel(combined_pb_de, combined_pb_null, "^PB_",
                              SCRNA_COLORS_V3, "A",
                              "Pseudobulk scRNA-seq")
p6b <- build_modality_panel(combined_prot_de, combined_prot_null, "^$",
                              PROT_COLORS_V3, "B",
                              "Proteomics")

# Panel C — spike-in F1 bar (4 datasets)
prot_de_ext <- copy(combined_prot_de)
prot_de_ext[, unified_F1 := fifelse(!is.na(F1) & F1 > 0, F1,
                             fifelse(!is.na(F1_fdr5) & F1_fdr5 > 0, F1_fdr5,
                                      NA_real_))]
spike_datasets <- c("sgsds_ratio2", "sgsds_ratio2.5", "ecoli_lfq", "cptac_ups1")
prot_spike <- prot_de_ext[dataset %in% spike_datasets & !is.na(unified_F1) & unified_F1 > 0]
prot_spike[, dataset_clean := clean_dataset_label(dataset)]
prot_spike[, method := factor(method, levels = names(PROT_COLORS_V3))]
prot_spike[, is_sgcb := method == "SGCB"]

p6c <- ggplot(prot_spike, aes(x = reorder(method, unified_F1, FUN = median),
                                 y = unified_F1, fill = method)) +
  geom_hline(yintercept = 0.5, color = "#D5DAE5", linetype = "dotted",
             linewidth = 0.26) +
  geom_hline(yintercept = 0.8, color = "#FF375F", linetype = "dashed",
             linewidth = 0.34) +
  geom_col(aes(alpha = is_sgcb), width = 0.72, color = "white",
           linewidth = 0.26) +
  geom_point(aes(size = is_sgcb), shape = 21, fill = "white",
             color = "#15151A", stroke = 0.3, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", unified_F1),
                fontface = ifelse(is_sgcb, "bold", "plain")),
             hjust = -0.12, size = 2.2, show.legend = FALSE) +
  facet_wrap(~dataset_clean, ncol = 4, scales = "free_x") +
  scale_fill_manual(values = PROT_COLORS_V3, drop = FALSE, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.58, `TRUE` = 1), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.45, `TRUE` = 2.35), guide = "none") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.5, 1)) +
  labs(title = "**C.**",
        x = NULL, y = "F1 @ FDR 5%") +
  THEME_JBHI(base_size = 10) +
  theme(axis.text.y = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         plot.margin = margin(2, 4, 4, 4))

fig6 <- cowplot::plot_grid(
  p6a, p6b, p6c,
  ncol = 1, rel_heights = c(0.98, 1.36, 0.78))
SAVE_FIG_V5("Fig6_multiomics", fig6, 7.5, 10.5)

# =========================================================================
# Fig 7 SEQC — facet lollipop polished
# =========================================================================
cat("  Fig 7 SEQC (v5)...\n")
seqc8 <- copy(combined_seqc)
seqc8[is.na(AUC) | (is.numeric(AUC) & is.nan(AUC)), AUC := 0]
seqc8[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
seqc_long <- melt(seqc8, id.vars = "method",
                   measure.vars = intersect(c("AUC", "F1", "precision"),
                                             names(seqc8)),
                   variable.name = "metric", value.name = "value")
seqc_long <- seqc_long[!(method == "samr")]
seqc_long[, metric := factor(metric, levels = c("AUC", "F1", "precision"),
                               labels = c("AUC", "F1", "Precision"))]
seqc_long[, is_best := value == max(value, na.rm = TRUE), by = metric]
seqc_long[, method := droplevels(method)]
seqc_long[, method := factor(method, levels = seqc8[!method %in% "samr"][order(-AUC)]$method)]

fig7_seqc <- ggplot(seqc_long, aes(x = value, y = method, color = method)) +
  geom_segment(aes(x = 0, xend = value, yend = method),
                linewidth = 0.35, alpha = 0.55) +
  geom_point(aes(size = is_best), show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", value),
                 fontface = ifelse(is_best, "bold", "plain")),
             hjust = -0.12, size = 2.4, show.legend = FALSE) +
  facet_wrap(~metric, ncol = 3, scales = "free_x") +
  scale_color_manual(values = METHOD_COLORS_V3, guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.8, `TRUE` = 4.5), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = NULL,
        x = "Value", y = NULL) +
  THEME_JBHI(base_size = 10)
SAVE_FIG_V5("Fig7_seqc_validation", fig7_seqc, 7.5, 4)

cat("  37b_v5: Benchmark figures done.\n")
