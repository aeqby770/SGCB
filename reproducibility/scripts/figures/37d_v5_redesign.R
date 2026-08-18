# =========================================================================
# 37d_v5_redesign.R — route-C figure redesign overlay
# 在 run_all 末尾 source（37b/37c 之后），用相同文件名覆盖被重设计的图。
# 直接复用 37b/37c 已加载到全局的 combined_* / kirc_* / cox_* 等数据。
# 设计基调（路线C）：克制骨架 + 关键处花哨 + SGCB 固定红高亮。
# =========================================================================
cat("\n=== 37d: route-C redesign overlay ===\n")

`%||%` <- function(a, b) if (is.null(a)) b else a
.md_sgcb <- function(levs)
  setNames(vapply(levs, function(m)
    if (m == "SGCB") "<span style='color:#FF3B30'>**SGCB**</span>" else m,
    character(1)), levs)

# 紧凑数据集短名：避免斜排长标签互相重叠 / 出界
.short_ds <- function(x) {
  m <- c(
    "cptac"           = "CPTAC",        "cptac_ups1"     = "CPTAC-UPS1",
    "ecoli_lfq"       = "E.coli LFQ",   "obrien_3species"= "O'Brien 3sp",
    "sgsds_ratio2"    = "SGSDS R2",     "sgsds_ratio2.5" = "SGSDS R2.5",
    "sim_proteomics"  = "Sim prot",     "tmt_mlr"        = "TMT MLR",
    "jakel"           = "Jakel MS",     "kang18"         = "Kang18",
    "muscat_sim"      = "Muscat sim",   "segerstolpe"    = "Segerstolpe",
    "xin"             = "Xin",          "zilionis"       = "Zilionis")
  out <- m[x]; out[is.na(out)] <- clean_dataset_label(x[is.na(out)]); unname(out)
}

# =========================================================================
# Fig 2 — Null calibration (A/B fingerprints kept; C de-cluttered)
# =========================================================================
redesign_fig2 <- function() {
  null_n10 <- copy(combined_null_pvalues)
  null_n10[, method := factor(method, levels = METHOD_ORDER_V3)]

  null_bins <- null_n10[, .N, by = .(method,
    bin = pmin(49L, as.integer(floor(pvalue * 50))))]
  null_bins[, total := sum(N), by = method]
  null_bins[, density := N / (total * 0.02)]
  null_bins[, bin_mid := (bin + 0.5) / 50]
  null_bins[, method := factor(method, levels = rev(METHOD_ORDER_V3))]
  ks_dt <- null_n10[, .(ks = round(ks.test(pvalue, "punif")$statistic, 3)),
                    by = method]
  ks_dt[, method := factor(method, levels = rev(METHOD_ORDER_V3))]
  ks_dt[, label := sprintf("KS %.3f", ks)]

  p2a <- ggplot(null_bins, aes(bin_mid, method, fill = density)) +
    geom_tile(width = 0.020, height = 0.78, color = "white", linewidth = 0.08) +
    geom_vline(xintercept = c(0.05, 0.95), color = "white", linewidth = 0.28,
               linetype = "dashed") +
    geom_text(data = ks_dt, aes(x = 1.04, y = method, label = label),
              inherit.aes = FALSE, hjust = 0, size = 2.25, color = "#1B1B1F",
              fontface = "bold") +
    scale_fill_gradientn(
      colors = c("#F8FBFF", "#C7F2EA", "#34D6C8", "#3157E8", "#281447"),
      limits = c(0, quantile(null_bins$density, 0.985, na.rm = TRUE)),
      oob = scales::squish, name = "Density") +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.18),
                       expand = expansion(0, 0)) +
    labs(title = "**A.** Null p-value fingerprint", x = "p-value", y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.y = ggtext::element_markdown(face = "bold", size = 8.4),
          legend.position = "bottom", legend.key.width = unit(1.05, "cm"),
          plot.margin = margin(4, 6, 4, 8)) +
    scale_y_discrete(labels = .md_sgcb(levels(null_bins$method)))

  qq_dev <- null_n10[, {
    n <- .N; theoretical <- ppoints(n); observed <- sort(pvalue)
    data.table(qbin = pmin(99L, as.integer(floor(theoretical * 100))),
               deviation = observed - theoretical)[, .(dev = mean(deviation)),
               by = qbin]
  }, by = method]
  qq_dev[, q_mid := (qbin + 0.5) / 100]
  qq_dev[, method := factor(method, levels = rev(METHOD_ORDER_V3))]
  p2b <- ggplot(qq_dev, aes(q_mid, method, fill = dev)) +
    geom_tile(width = 0.010, height = 0.78, color = NA) +
    geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "white",
               linewidth = 0.18, alpha = 0.75) +
    scale_fill_gradient2(low = "#E9435D", mid = "#F9FAFC", high = "#2156D9",
      midpoint = 0, limits = c(-0.35, 0.35), oob = scales::squish,
      name = "Observed - expected") +
    scale_x_continuous(breaks = c(0, 0.5, 1), expand = expansion(0, 0)) +
    labs(title = "**B.** Q-Q deviation field",
         x = "Theoretical quantile", y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.y = ggtext::element_markdown(face = "bold", size = 8.4),
          legend.position = "bottom", legend.key.width = unit(1.05, "cm"),
          plot.margin = margin(4, 6, 4, 8)) +
    scale_y_discrete(labels = .md_sgcb(levels(qq_dev$method)))

  # Panel C — de-cluttered: thin error bars (no overlapping ribbons)
  fdr_cal8 <- melt(combined_sim[, .(method, actual_FDR_fdr1, actual_FDR_fdr5,
                                     actual_FDR_fdr10)],
                   id.vars = "method",
                   measure.vars = c("actual_FDR_fdr1", "actual_FDR_fdr5",
                                    "actual_FDR_fdr10"),
                   variable.name = "threshold", value.name = "observed_fdr")
  fdr_cal8[, nominal := fifelse(threshold == "actual_FDR_fdr1", 0.01,
                         fifelse(threshold == "actual_FDR_fdr5", 0.05, 0.10))]
  fdr_summary <- fdr_cal8[, .(mean_fdr = mean(observed_fdr, na.rm = TRUE),
                              se_fdr = sd(observed_fdr, na.rm = TRUE) /
                                       sqrt(.N)),
                          by = .(method, nominal)]
  fdr_summary[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
  fdr_summary[is.na(se_fdr), se_fdr := 0]
  calibrated <- c("SGCB", "limma", "glmGamPoi")
  fdr_summary[, group := fifelse(method %in% calibrated, "Calibrated",
                                 "Anti-conservative")]
  fdr_summary[, is_sgcb := method == "SGCB"]
  p2c <- ggplot(fdr_summary, aes(nominal, mean_fdr, color = method)) +
    annotate("rect", xmin = 0, xmax = 0.105, ymin = 0, ymax = 0.05,
             fill = "#37C98E", alpha = 0.10) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45",
                linewidth = 0.35) +
    geom_line(aes(linewidth = is_sgcb), alpha = 0.92) +
    geom_pointrange(aes(ymin = pmax(mean_fdr - se_fdr, 0),
                        ymax = mean_fdr + se_fdr, size = is_sgcb),
                    fatten = 1.9, linewidth = 0.45, shape = 21, fill = "white",
                    stroke = 0.8) +
    facet_wrap(~group, scales = "free_y") +
    scale_color_manual(values = METHOD_COLORS_V3, name = "Method") +
    scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.0),
                           guide = "none") +
    scale_size_manual(values = c(`FALSE` = 0.28, `TRUE` = 0.5), guide = "none") +
    scale_x_continuous(breaks = c(0.01, 0.05, 0.10),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "**C.** Calibration operating zone",
         x = "Nominal FDR", y = "Observed FDR") +
    THEME_JBHI(base_size = 10) +
    theme(legend.key.size = unit(0.35, "cm"),
          strip.text = element_text(face = "bold"),
          plot.margin = margin(4, 8, 4, 8))

  fig2 <- cowplot::plot_grid(p2a, p2b, p2c, ncol = 1,
                             rel_heights = c(0.86, 0.86, 1.0))
  SAVE_FIG_V5("Fig2_calibration", fig2, 7.5, 8.2)
}

# =========================================================================
# Fig 3 — Simulation power (heatmap WITHOUT side-bar; FDR dumbbell-ish)
# =========================================================================
redesign_fig3 <- function() {
  sim_agg <- combined_sim[, .(mean_F1 = F1_fdr5, mean_FDR = actual_FDR_fdr5,
                              n_per_group = as.integer(sub("sim_n([0-9]+)_.*",
                                                           "\\1", dataset)),
                              effect_size = as.numeric(sub(".*_ef([0-9.]+)$",
                                                            "\\1", dataset))),
                          by = .(method, dataset)]
  sim_agg[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
  sim_agg[, ef_lab := paste0("Effect = ", effect_size, "\u00d7")]
  ord_tbl <- sim_agg[n_per_group == 10 & effect_size == 2,
                     .(method, mean_F1)][order(-mean_F1)]
  meth_ord <- as.character(ord_tbl$method)
  sim_agg[, method := factor(method, levels = rev(meth_ord))]
  sim_agg[, F1_txt_col := ifelse(mean_F1 > 0.45, "white", "grey15")]
  sgcb_tiles <- sim_agg[as.character(method) == "SGCB"]

  p3a <- ggplot(sim_agg, aes(factor(n_per_group), method, fill = mean_F1)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_tile(data = sgcb_tiles, aes(factor(n_per_group), method),
              inherit.aes = FALSE, fill = NA, color = "#FF3B30",
              linewidth = 0.75) +
    geom_text(aes(label = sprintf("%.2f", mean_F1), color = F1_txt_col),
              size = 2.7, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(colors = c("#FAFCFF", "#D9FFF7", "#52E4CF",
                                    "#2A70FF", "#211442"),
                         limits = c(0, 1), name = "F1", breaks = c(0, 0.5, 1),
                         guide = guide_colorbar(barwidth = 10, barheight = 0.4,
                           title.position = "left", title.vjust = 1)) +
    facet_wrap(~ef_lab, ncol = 2) +
    scale_x_discrete(expand = expansion(0, 0)) +
    scale_y_discrete(expand = expansion(0, 0),
                     labels = .md_sgcb(levels(sim_agg$method))) +
    labs(title = "**A.** F1 across sample size & effect",
         x = "Samples per group", y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(strip.text = element_text(size = 8.5, face = "bold"),
          legend.position = "bottom", legend.justification = "left",
          panel.border = element_rect(color = "grey30", fill = NA,
                                      linewidth = 0.3),
          axis.text.y = ggtext::element_markdown(face = "bold", size = 8.5))

  # Panel B — FDR control as dumbbell (1.5x -> 2x) per method
  fdr_db <- dcast(sim_agg[, .(fdr = mean(mean_FDR, na.rm = TRUE)),
                          by = .(method, effect_size, n_per_group)],
                  method + n_per_group ~ effect_size, value.var = "fdr")
  setnames(fdr_db, c("1.5", "2"), c("e15", "e20"), skip_absent = TRUE)
  fdr_db[, method := factor(as.character(method), levels = rev(meth_ord))]
  fdr_db[, n_lab := factor(paste0("n=", n_per_group),
                           levels = paste0("n=", sort(unique(n_per_group))))]
  fdr_db_long <- melt(fdr_db, id.vars = c("method", "n_per_group", "n_lab"),
                      measure.vars = c("e15", "e20"),
                      variable.name = "effect", value.name = "fdr")
  fdr_db_long[, fdr := pmin(fdr, 0.25)]
  fdr_db_long[, effect := factor(effect, levels = c("e15", "e20"),
                                 labels = c("1.5x", "2x"))]
  p3b <- ggplot(fdr_db, aes(y = method)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "#FF1744",
               linewidth = 0.35) +
    geom_segment(aes(x = pmin(e15, 0.25), xend = pmin(e20, 0.25),
                     yend = method), color = "grey70", linewidth = 0.9) +
    geom_point(data = fdr_db_long, aes(x = fdr, shape = effect, color = method),
               size = 2.4) +
    facet_wrap(~n_lab, ncol = 1) +
    scale_color_manual(values = METHOD_COLORS_V3, guide = "none") +
    scale_shape_manual(values = c(`1.5x` = 16, `2x` = 17),
                       name = "Effect") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 0.26), expand = expansion(0, 0)) +
    labs(title = "**B.** Observed FDR (1.5\u00d7 \u2192 2\u00d7)",
         x = "Observed FDR", y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(legend.position = "bottom",
          axis.text.y = ggtext::element_markdown(size = 7.5),
          strip.text = element_text(face = "bold", size = 8))

  fig3 <- cowplot::plot_grid(p3a, p3b, ncol = 1, rel_heights = c(1.05, 1.3))
  SAVE_FIG_V5("Fig3_simulation_power", fig3, 7.5, 8.6)
}

# =========================================================================
# Fig 4 — Practical (A yield heatmap + Mean column; B SE pointrange; C FDR)
# =========================================================================
redesign_fig4 <- function() {
  classic8 <- copy(combined_classic)
  classic8[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
  cs <- classic8[, .(n_sig = sum(n_sig_fdr5, na.rm = TRUE)),
                 by = .(method, dataset)]
  cs[, dataset_clean := clean_dataset_label(dataset)]
  cmean <- cs[, .(dataset_clean = "Mean", n_sig = mean(n_sig, na.rm = TRUE)),
              by = method]
  cs2 <- rbind(cs[, .(method, dataset_clean, n_sig)], cmean)
  cs2[, log_n := fifelse(n_sig <= 0, 0, log10(n_sig))]
  ds_levels <- c(sort(unique(cs$dataset_clean)), "Mean")
  cs2[, dataset_clean := factor(dataset_clean, levels = ds_levels)]
  ord <- cmean[order(n_sig)]$method
  cs2[, method := factor(as.character(method), levels = ord)]
  cs2[, txt_col := ifelse(log_n > 2.5, "white", "grey15")]
  sgcb_t <- cs2[as.character(method) == "SGCB"]

  p4a <- ggplot(cs2, aes(dataset_clean, method, fill = log_n)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_tile(data = cs2[dataset_clean == "Mean"],
              aes(dataset_clean, method), inherit.aes = FALSE, fill = NA,
              color = "grey55", linewidth = 0.4) +
    geom_tile(data = sgcb_t, aes(dataset_clean, method), inherit.aes = FALSE,
              fill = NA, color = "#FF3B30", linewidth = 0.7) +
    geom_text(aes(label = ifelse(n_sig == 0, "0", scales::comma(round(n_sig))),
                  color = txt_col), size = 2.4, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(colors = c("#FAFCFF", "#D9FFF7", "#52E4CF",
                                    "#2A70FF", "#211442"),
                         limits = c(0, max(cs2$log_n) * 1.01),
                         name = expression(log[10]*"(DE genes)"),
                         guide = guide_colorbar(barwidth = 0.4, barheight = 6,
                                                title.position = "top")) +
    scale_x_discrete(expand = expansion(0, 0)) +
    scale_y_discrete(expand = expansion(0, 0),
                     labels = .md_sgcb(levels(cs2$method))) +
    labs(title = "**A.** DE yield on classic datasets", x = NULL, y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          axis.text.y = ggtext::element_markdown(face = "bold", size = 8.5),
          legend.position = "right",
          panel.border = element_rect(color = "grey30", fill = NA,
                                      linewidth = 0.3))

  # Panel B — F1 across n with SE pointrange (replaces glud SD ribbon)
  ss <- copy(combined_sim)
  ss[, n_per_group := as.integer(sub("sim_n([0-9]+)_.*", "\\1", dataset))]
  ss <- ss[, .(mean_F1 = mean(F1_fdr5, na.rm = TRUE),
               se_F1 = sd(F1_fdr5, na.rm = TRUE) / sqrt(.N)),
           by = .(method, n_per_group)]
  ss[is.na(se_F1), se_F1 := 0]
  ss[, method := factor(as.character(method), levels = METHOD_ORDER_V3)]
  ss[, is_sgcb := method == "SGCB"]
  p4b <- ggplot(ss, aes(factor(n_per_group), mean_F1, color = method,
                        group = method)) +
    geom_line(aes(linewidth = is_sgcb, alpha = is_sgcb),
              position = position_dodge(width = 0.4)) +
    geom_pointrange(aes(ymin = pmax(mean_F1 - se_F1, 0),
                        ymax = pmin(mean_F1 + se_F1, 1), size = is_sgcb,
                        alpha = is_sgcb),
                    position = position_dodge(width = 0.4),
                    linewidth = 0.5, shape = 21, fill = "white",
                    stroke = 0.7) +
    scale_color_manual(values = METHOD_COLORS_V3, name = "Method") +
    scale_alpha_manual(values = c(`FALSE` = 0.38, `TRUE` = 1), guide = "none") +
    scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.1),
                           guide = "none") +
    scale_size_manual(values = c(`FALSE` = 0.22, `TRUE` = 0.45), guide = "none") +
    labs(title = "**B.** Mean F1 \u00b1 SE across sample size",
         x = "Samples per group", y = "Mean F1") +
    THEME_JBHI(base_size = 10) + theme(legend.position = "right")

  # Panel C — FDR heatmap (kept)
  ssc <- copy(ss); ssc[, mean_FDR := NA_real_]
  ssc2 <- copy(combined_sim)
  ssc2[, n_per_group := as.integer(sub("sim_n([0-9]+)_.*", "\\1", dataset))]
  ssc2 <- ssc2[, .(mean_FDR = mean(actual_FDR_fdr5, na.rm = TRUE)),
               by = .(method, n_per_group)]
  ssc2[, method := factor(as.character(method), levels = rev(ord))]
  ssc2[, txt_col := ifelse(mean_FDR > 0.08, "white", "grey15")]
  p4c <- ggplot(ssc2, aes(factor(n_per_group), method, fill = mean_FDR)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", mean_FDR), color = txt_col),
              size = 2.4, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(colors = c("#FAFCFF", "#B9F6E3", "#FFCA7A",
                                    "#FF5A5F", "#5B1022"),
                         values = scales::rescale(c(0, 0.05, 0.2, 0.7)),
                         limits = c(0, max(ssc2$mean_FDR, na.rm = TRUE) * 1.01),
                         name = "Obs. FDR",
                         guide = guide_colorbar(barwidth = 0.4, barheight = 5,
                                                title.position = "top")) +
    scale_y_discrete(labels = .md_sgcb(levels(ssc2$method))) +
    labs(title = "**C.** Observed FDR across sample size",
         x = "Samples per group", y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.y = ggtext::element_markdown(face = "bold", size = 8),
          legend.position = "right",
          panel.border = element_rect(color = "grey30", fill = NA,
                                      linewidth = 0.3))

  p4bc <- cowplot::plot_grid(p4b, p4c, ncol = 2, rel_widths = c(1.4, 1),
                             align = "h", axis = "tb")
  fig4 <- cowplot::plot_grid(p4a, p4bc, ncol = 1, rel_heights = c(1.25, 1))
  SAVE_FIG_V5("Fig4_practical", fig4, 7.5, 9.3)
}

# =========================================================================
# Fig 6 — Multiomics bubble matrix (size=DE genes, fill=null FDR)
# =========================================================================
.bubble_panel <- function(de_dt, null_dt, strip_re, letter, label) {
  de <- copy(de_dt)[!grepl("^dv_", method) & !is.na(method) & method != ""]
  nu <- copy(null_dt)[!grepl("^dv_", method) & !is.na(method) & method != ""]
  de[, ds := .short_ds(gsub(strip_re, "", dataset))]
  de <- de[, .(n_sig = sum(n_sig_fdr5, na.rm = TRUE)), by = .(method, ds)]
  de[, log_n := fifelse(n_sig <= 0, 0, log10(n_sig))]
  num <- nu[, .(fdr = mean(mean_FDR05, na.rm = TRUE)), by = method]
  de <- merge(de, num, by = "method", all.x = TRUE)
  ord <- de[, .(m = mean(n_sig, na.rm = TRUE)), by = method][order(m)]$method
  de[, method := factor(method, levels = ord)]
  de[, over := !is.na(fdr) & fdr > 0.05]
  maxf <- max(de$fdr, 0.06, na.rm = TRUE)
  zebra <- data.table(yi = seq_along(levels(de$method)))[yi %% 2 == 0]

  ggplot(de, aes(ds, method)) +
    geom_rect(data = zebra, aes(ymin = yi - 0.5, ymax = yi + 0.5),
              xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#EEF2F8") +
    geom_point(aes(size = log_n, fill = fdr), shape = 21, color = "grey25",
               stroke = 0.34) +
    geom_point(data = de[over == TRUE], aes(size = log_n), shape = 21,
               fill = NA, color = "#B2182B", stroke = 1.0) +
    scale_size_continuous(range = c(0.7, 8.4), breaks = c(1, 2, 3, 4),
      labels = c("10", "100", "1k", "10k"), name = "DE genes",
      limits = c(0, NA)) +
    scale_fill_gradientn(
      colors = c("#2166AC", "#7FB1D6", "#F7F7F7", "#EF8A62", "#B2182B"),
      values = scales::rescale(c(0, 0.025, 0.05, 0.5 * maxf, maxf),
                               from = c(0, maxf)),
      limits = c(0, maxf), oob = scales::squish,
      breaks = c(0, 0.05, round(maxf, 2)), name = "Null FDR",
      guide = guide_colorbar(barwidth = 0.5, barheight = 4.8,
                             title.position = "top", order = 2,
                             frame.colour = "grey40", ticks.colour = "grey40")) +
    scale_y_discrete(labels = .md_sgcb(levels(de$method)),
                     expand = expansion(0, 0.6)) +
    scale_x_discrete(position = "top", expand = expansion(0, 0.62)) +
    labs(title = sprintf("**%s.** %s", letter, label), x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0,
                                         size = 8.0, color = "grey15"),
          axis.text.y = ggtext::element_markdown(size = 8.6),
          plot.title = ggtext::element_markdown(size = 11, face = "bold",
                                                margin = margin(b = 14)),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.26),
          panel.border = element_rect(color = "grey35", fill = NA,
                                      linewidth = 0.35),
          legend.position = "right", legend.box = "vertical",
          legend.spacing.y = unit(2, "pt"),
          plot.margin = margin(10, 10, 4, 6))
}

redesign_fig6 <- function() {
  p6a <- .bubble_panel(combined_pb_de, combined_pb_null, "^PB_", "A",
                       "Pseudobulk scRNA-seq")
  p6b <- .bubble_panel(combined_prot_de, combined_prot_null, "^$", "B",
                       "Proteomics")
  # Panel C — spike-in F1: lollipop + 数字右对齐成列(防重叠/出界) + best 描圈
  prot <- copy(combined_prot_de)
  prot[, uF1 := fifelse(!is.na(F1) & F1 > 0, F1,
                fifelse(!is.na(F1_fdr5) & F1_fdr5 > 0, F1_fdr5, NA_real_))]
  spk <- c("cptac_ups1", "ecoli_lfq", "sgsds_ratio2", "sgsds_ratio2.5")
  prot <- prot[dataset %in% spk & !is.na(uF1) & uF1 > 0]
  prot[, ds := factor(.short_ds(dataset), levels = .short_ds(spk))]
  ord <- prot[, .(med = median(uF1, na.rm = TRUE)), by = method][order(med)]$method
  prot[, method := factor(method, levels = ord)]
  prot[, is_sgcb := method == "SGCB"]
  prot[, is_best := uF1 == max(uF1), by = ds]
  ylv <- levels(prot$method)
  zebra <- data.table(yi = seq_along(ylv))[yi %% 2 == 0]
  p6c <- ggplot(prot, aes(x = uF1, y = method)) +
    geom_rect(data = zebra, aes(ymin = yi - 0.5, ymax = yi + 0.5),
              xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#EEF2F8") +
    geom_vline(xintercept = 0.8, color = "#FF375F", linetype = "dashed",
               linewidth = 0.36) +
    geom_segment(aes(x = 0, xend = uF1, yend = method, color = method),
                 linewidth = 0.7, alpha = 0.55) +
    geom_point(aes(color = method, size = is_sgcb)) +
    geom_point(data = prot[is_best == TRUE], shape = 21, fill = NA,
               color = "#1B1B1F", stroke = 0.7, size = 3.6) +
    geom_text(data = prot[is_sgcb == FALSE],
              aes(x = Inf, label = sprintf("%.2f", uF1)),
              hjust = 1.0, size = 2.6, color = "#1B1B1F") +
    geom_text(data = prot[is_sgcb == TRUE],
              aes(x = Inf, label = sprintf("%.2f", uF1)),
              hjust = 1.0, size = 2.7, color = "#FF3B30", fontface = "bold") +
    facet_wrap(~ds, ncol = 4) +
    scale_color_manual(values = PROT_COLORS_V3, guide = "none") +
    scale_size_manual(values = c(`FALSE` = 2.1, `TRUE` = 3.6), guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0.02, 0.22))) +
    scale_y_discrete(labels = .md_sgcb(ylv)) +
    labs(title = "**C.** Spike-in F1 @ FDR 5% (ground-truth datasets)",
         x = "F1", y = NULL) +
    coord_cartesian(clip = "off") +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.y = ggtext::element_markdown(size = 8.0),
          plot.title = ggtext::element_markdown(size = 11, face = "bold",
                                                margin = margin(b = 6)),
          strip.text = element_text(size = 8.4, face = "bold"),
          panel.spacing.x = unit(0.85, "lines"),
          panel.grid = element_blank(),
          plot.margin = margin(6, 12, 4, 6))

  fig6 <- cowplot::plot_grid(p6a, p6b, p6c, ncol = 1,
                             rel_heights = c(1.32, 0.94, 0.78))
  SAVE_FIG_V5("Fig6_multiomics", fig6, 7.6, 11.0)
}

# =========================================================================
# Fig 7 SEQC — parallel coordinates (replaces triple lollipop)
# =========================================================================
redesign_fig7_seqc <- function() {
  seqc <- copy(combined_seqc)
  seqc[is.na(AUC) | is.nan(AUC), AUC := 0]
  seqc <- seqc[method != "samr"]
  mets <- c("AUC", "F1", "precision")
  long <- melt(seqc, id.vars = "method", measure.vars = mets,
               variable.name = "metric", value.name = "value")
  long[, metric := factor(metric, levels = mets,
                          labels = c("AUC", "F1", "Precision"))]
  long[, vnorm := (value - min(value)) / (max(value) - min(value)), by = metric]
  long[, is_sgcb := method == "SGCB"]
  ord <- seqc[order(-AUC)]$method
  long[, method := factor(method, levels = ord)]
  axis_ann <- long[, .(vmin = min(value), vmax = max(value)), by = metric]
  axis_ann[, xi := as.integer(metric)]
  long[, is_best := value == max(value), by = metric]

  p <- ggplot(long, aes(metric, vnorm, group = method)) +
    geom_vline(xintercept = 1:3, color = "grey80", linewidth = 0.4) +
    geom_line(data = long[is_sgcb == FALSE], aes(color = method),
              linewidth = 0.55, alpha = 0.65) +
    geom_point(data = long[is_sgcb == FALSE], aes(color = method),
               size = 1.9, alpha = 0.8) +
    geom_line(data = long[is_sgcb == TRUE], color = "#FF3B30", linewidth = 1.7) +
    geom_point(data = long[is_sgcb == TRUE], color = "#FF3B30", size = 3.4) +
    geom_point(data = long[is_best == TRUE], shape = 21, fill = NA,
               color = "grey15", size = 4.4, stroke = 0.6) +
    ggrepel::geom_text_repel(data = long[metric == "Precision"],
      aes(label = sprintf("%s  %.3f", method, value), color = method),
      hjust = 0, direction = "y", nudge_x = 0.20, segment.size = 0.2,
      size = 2.5, show.legend = FALSE, max.overlaps = Inf,
      box.padding = 0.55, point.padding = 0.2, force = 8,
      min.segment.length = 0, seed = 1) +
    ggrepel::geom_text_repel(data = long[metric == "AUC"],
      aes(label = sprintf("%.3f", value), color = method),
      hjust = 1, direction = "y", nudge_x = -0.16, segment.size = 0.2,
      size = 2.3, show.legend = FALSE, max.overlaps = Inf,
      box.padding = 0.5, point.padding = 0.2, force = 8,
      min.segment.length = 0, seed = 1) +
    geom_text(data = axis_ann, aes(x = xi, y = 1.07, label = sprintf("%.2f", vmax)),
              inherit.aes = FALSE, size = 2.3, color = "grey35") +
    geom_text(data = axis_ann, aes(x = xi, y = -0.07, label = sprintf("%.2f", vmin)),
              inherit.aes = FALSE, size = 2.3, color = "grey35") +
    scale_color_manual(values = METHOD_COLORS_V3, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.30, 0.42))) +
    scale_y_continuous(limits = c(-0.12, 1.14), expand = expansion(0, 0)) +
    labs(title = "**SEQC / MAQC-III external validation**",
         subtitle = "Each line = one method across three ground-truth metrics; SGCB highlighted",
         x = NULL, y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(face = "bold", size = 10),
          plot.title = ggtext::element_markdown(size = 11, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 8, color = "grey40"),
          plot.margin = margin(8, 10, 6, 10))
  SAVE_FIG_V5("Fig7_seqc_validation", p, 7.5, 4.2)
}

.has <- function(...) all(vapply(c(...), exists, logical(1), envir = .GlobalEnv))

# =========================================================================
# Fig 7 — TCGA-KIRC merged: re-assemble to use p7e (fills KM whitespace)
# =========================================================================
redesign_fig7 <- function() {
  if (!.has("p7a", "p7b", "p7c", "p7d", "p7e")) {
    cat("  [37d] Fig7 panels missing; skip\n"); return(invisible()) }
  # rebuild Panel E with short, non-truncated facet titles
  p7e_use <- p7e
  if (exists("expr_long", envir = .GlobalEnv)) {
    el <- copy(expr_long)
    el[, gshort := sub(" \\(.*$", "", as.character(facet_lbl))]
    el[, gshort := factor(gshort, levels = unique(gshort[order(gene_order)]))]
    p7e_use <- ggplot(el, aes(x = group, y = expr, fill = group)) +
      ggdist::stat_halfeye(adjust = 0.8, width = 0.5, .width = 0,
                           justification = -0.15, point_colour = NA, alpha = 0.5) +
      geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 0.9,
                   color = "grey25", linewidth = 0.3) +
      facet_wrap(~gshort, ncol = 6, scales = "free_y") +
      scale_fill_manual(values = c(normal = "#4DAF4A", tumor = "#E41A1C"),
                        guide = "none") +
      labs(title = "**E.** Representative channel genes (DE / DV)",
           x = NULL, y = expression(log[2] ~ "expression")) +
      THEME_JBHI(base_size = 10) +
      theme(strip.text = element_text(size = 7.6, face = "bold"),
            axis.text.x = element_text(size = 7),
            panel.spacing = unit(0.3, "lines"),
            plot.title = ggtext::element_markdown(size = 10.5, face = "bold"))
  }
  row1 <- cowplot::plot_grid(p7a, p7b, p7c, ncol = 3,
                             rel_widths = c(1, 1.2, 1.3), align = "h", axis = "tb")
  fig7 <- cowplot::plot_grid(row1, p7d, p7e_use, ncol = 1,
                             rel_heights = c(1.0, 0.82, 0.62))
  SAVE_FIG_V5("Fig7_tcga_kirc_merged", fig7, 8.0, 9.4)
}

# =========================================================================
# Fig 8 — kang18: degenerate donut (99% DE-only) -> log-scale bar
# =========================================================================
redesign_fig8 <- function() {
  if (!exists("kang_ch_plot", envir = .GlobalEnv) ||
      !.has("p8b", "p8c", "p8d")) {
    cat("  [37d] Fig8 inputs missing; skip\n"); return(invisible()) }
  d <- copy(kang_ch_plot)
  KCOL <- c("DE-only" = "#2F5BFF", "DV-only" = "#15C8A2",
            "DG-only" = "#8A4DFF", "DE & DV" = "#FFB000")
  d[, cs := factor(channel_short, levels = c("DE-only", "DV-only",
                                             "DG-only", "DE & DV"))]
  p8a_new <- ggplot(d, aes(x = n_features, y = reorder(cs, n_features),
                           fill = cs)) +
    geom_col(width = 0.66, color = "white", linewidth = 0.4) +
    geom_text(aes(label = scales::comma(n_features)), hjust = -0.2,
              size = 2.8, fontface = "bold", color = "#15151A") +
    scale_fill_manual(values = KCOL, guide = "none") +
    scale_x_log10(expand = expansion(mult = c(0, 0.32)),
                  labels = scales::comma,
                  breaks = c(1, 10, 100, 1000)) +
    annotate("text", x = 1, y = 0.4,
             label = sprintf("total = %s active genes",
                             scales::comma(sum(d$n_features))),
             hjust = 0, vjust = 0, size = 2.4, color = "grey45") +
    labs(title = "**A.** Channel composition (log scale)",
         x = "Genes", y = NULL) +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.y = element_text(face = "bold", size = 8.5),
          plot.title = ggtext::element_markdown(size = 10.5, face = "bold"),
          panel.grid.major.y = element_blank())
  # rebuild Panel C so the n= labels are not covered by large points
  p8c_use <- p8c
  if (exists("kang_mk2", envir = .GlobalEnv)) {
    p8c_use <- ggplot(kang_mk2, aes(x = fold_enrichment, y = cell_type_wrap)) +
      geom_segment(aes(x = 0, xend = fold_enrichment, yend = cell_type_wrap),
                   linewidth = 0.4, alpha = 0.6, color = "#3B4CC0") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey55",
                 linewidth = 0.3) +
      geom_point(aes(size = neglog10p), shape = 21, fill = "#3B4CC0",
                 color = "grey20", stroke = 0.3, alpha = 0.92) +
      geom_text(aes(label = n_lbl), hjust = 0, nudge_x = 0.25, size = 2.2,
                color = "grey25") +
      scale_size_continuous(range = c(2.4, 6.0),
                            name = expression(-log[10] ~ q)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.38))) +
      labs(title = "**C.**", x = "Fold enrichment", y = NULL) +
      THEME_JBHI(base_size = 10) +
      theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
            axis.text.y = element_text(size = 7),
            plot.title = ggtext::element_markdown(size = 10.5, face = "bold"),
            plot.margin = margin(2, 6, 4, 4))
  }
  top_8 <- cowplot::plot_grid(p8a_new, p8b, ncol = 2,
                              rel_widths = c(0.95, 1.3), align = "h", axis = "tb")
  bottom_8 <- cowplot::plot_grid(p8c_use, p8d, ncol = 2,
                                 rel_widths = c(0.82, 1.43), align = "h", axis = "tb")
  fig8 <- cowplot::plot_grid(top_8, bottom_8, ncol = 1,
                             rel_heights = c(0.92, 1.04))
  SAVE_FIG_V5("Fig8_kang18_merged", fig8, 7.5, 8.7)
}

# =========================================================================
# FigS1 — manhattan: drop empty channels (DG-only n=0) to reclaim space
# =========================================================================
redesign_figS1 <- function() {
  if (!exists("kirc_dt_all", envir = .GlobalEnv) ||
      !exists("CH_COLS_V5", envir = .GlobalEnv)) {
    cat("  [37d] FigS1 inputs missing; skip\n"); return(invisible()) }
  kd <- copy(kirc_dt_all)
  cnts <- kd[, .N, by = channel]
  keep <- as.character(cnts[N > 0L]$channel)          # drop only truly-empty (DG-only n=0)
  keep <- intersect(names(CH_COLS_V5), keep)
  big_ch <- as.character(cnts[N >= 10L]$channel)      # violins only where enough points
  kd <- kd[as.character(channel) %in% keep]
  kd[, channel := factor(as.character(channel), levels = keep)]
  kd[, ch_num := as.integer(channel)]
  set.seed(12345)
  kd[, x_jit := ch_num + runif(.N, -0.36, 0.36)]
  y_cap <- quantile(abs(kd$log2FC), 0.995, na.rm = TRUE)
  kd[, y_show := pmin(pmax(log2FC, -y_cap), y_cap)]
  kd[, dist_score := abs(log2FC) * (-log10(pmax(padj_de, 1e-300)))]
  nlbl <- c("DE-only" = 5L, "DV-only" = 4L, "DE & DV" = 4L,
            "DD-only (not DE)" = 5L)
  top_genes <- unlist(lapply(keep, function(ch) {
    sub <- kd[channel == ch]; k <- nlbl[ch]; if (is.na(k)) k <- 4L
    if (nrow(sub) == 0L) return(character())
    head(sub[order(-dist_score)], k)$gene }))
  top_genes <- unique(top_genes[!is.na(top_genes)])
  CH_BG <- colorspace::lighten(CH_COLS_V5[keep], 0.56)
  bg_df <- data.table(xmin = seq_along(keep) - 0.5, xmax = seq_along(keep) + 0.5,
                      channel = factor(keep, levels = keep))
  cc <- merge(data.table(channel = factor(keep, levels = keep)),
              kd[, .N, by = channel], by = "channel", all.x = TRUE)
  cc[is.na(N), N := 0L]
  cc[, label := paste0(channel, "\n(n = ", scales::comma(N), ")")]
  lfc_lim <- max(abs(kd$y_show), na.rm = TRUE)
  figS1 <- ggplot(kd, aes(x = x_jit, y = y_show)) +
    geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax, ymin = -y_cap * 1.05,
              ymax = y_cap * 1.05, fill = channel), inherit.aes = FALSE,
              show.legend = FALSE) +
    geom_violin(data = kd[as.character(channel) %in% big_ch],
                aes(x = ch_num, y = y_show, group = channel,
                fill = channel), inherit.aes = FALSE, width = 0.88, alpha = 0.42,
                color = "white", linewidth = 0.36, trim = TRUE) +
    geom_hline(yintercept = 0, color = "grey25", linewidth = 0.3) +
    ggrastr::rasterise(geom_point(aes(color = y_show), size = 0.24, alpha = 0.34,
                       shape = 16), dpi = 300) +
    geom_point(data = kd[gene %in% top_genes], aes(color = y_show), size = 1.8,
               alpha = 0.95, shape = 16) +
    geom_point(data = kd[gene %in% top_genes], color = "black", size = 2.2,
               alpha = 0.8, shape = 1, stroke = 0.55) +
    ggrepel::geom_text_repel(data = kd[gene %in% top_genes], aes(label = gene),
      size = 2.5, color = "black", bg.color = "white", bg.r = 0.16,
      max.overlaps = Inf, segment.size = 0.22, segment.color = "grey30",
      box.padding = 0.6, point.padding = 0.35, force = 32, force_pull = 0.1,
      seed = 42, min.segment.length = 0, max.iter = 30000) +
    scale_color_gradientn(name = expression(log[2]*"FC"),
      colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
      limits = c(-lfc_lim, lfc_lim), oob = scales::squish,
      breaks = c(-lfc_lim, 0, lfc_lim), labels = function(x) sprintf("%.1f", x),
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
        barwidth = unit(0.5, "cm"), barheight = unit(3.2, "cm"),
        title.hjust = 0.5, title.position = "top")) +
    scale_fill_manual(values = CH_BG) +
    scale_x_continuous(breaks = seq_along(keep), labels = cc$label,
                       expand = expansion(mult = 0.02)) +
    coord_cartesian(ylim = c(-y_cap * 1.05, y_cap * 1.05), clip = "on") +
    labs(title = NULL, x = NULL,
         y = expression(Average ~ log[2] ~ "fold change")) +
    THEME_JBHI(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold", size = 9, color = "black",
                                     margin = margin(t = 5), lineheight = 1.05),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(color = "#15151A", fill = NA,
                                      linewidth = 0.35),
          legend.position = "right")
  SAVE_FIG_V5("FigS1_tcga_kirc_manhattan", figS1, 11.5, 6.5)
}

# =========================================================================
# FigS3 — pan-cancer: replace focus lollipop (B) with grouped bar
# =========================================================================
redesign_figS3 <- function() {
  if (!exists("pan_l", envir = .GlobalEnv) || !exists("pS3a", envir = .GlobalEnv) ||
      !exists("CH_COLS_V5", envir = .GlobalEnv)) {
    cat("  [37d] FigS3 inputs missing; skip\n"); return(invisible()) }
  pan_foc <- pan_l[channel %in% c("DV-only", "DD-only (not DE)")]
  pS3b_new <- ggplot(pan_foc, aes(x = cancer, y = count, fill = channel)) +
    geom_col(position = position_dodge(width = 0.74), width = 0.7,
             color = "white", linewidth = 0.3) +
    geom_text(aes(label = scales::comma(count)),
              position = position_dodge(width = 0.74), vjust = -0.4,
              size = 2.0, color = "grey25") +
    scale_fill_manual(values = CH_COLS_V5[c("DV-only", "DD-only (not DE)")],
                      name = "Channel", guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(title = "**B.** Variability-driven genes per cancer (<span style='color:#00A77E'>DV-only</span> / <span style='color:#FF6D00'>DD-only</span>)",
         x = NULL, y = "Genes") +
    THEME_JBHI(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = ggtext::element_markdown(size = 10))
  figS3 <- (pS3a / pS3b_new) + plot_layout(heights = c(1.2, 1),
            guides = "collect") &
    theme(legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
          legend.text = element_text(size = 7))
  SAVE_FIG_V5("FigS3_pan_cancer_channels", figS3, 10, 9)
}

# =========================================================================
# FigS5 — clinical: DEDUP vs Fig7 (drop forest+KM); show by-stage signatures
# =========================================================================
redesign_figS5 <- function() {
  if (!exists("matched", envir = .GlobalEnv) ||
      !exists("stage_tab", envir = .GlobalEnv)) {
    cat("  [37d] FigS5 inputs missing; skip\n"); return(invisible()) }
  sc <- grep("^score_", names(matched), value = TRUE)
  sc <- sc[vapply(sc, function(c) sum(!is.na(matched[[c]])) > 10L, logical(1))]
  if (length(sc) == 0L) { cat("  [37d] FigS5 no scores; skip\n"); return(invisible()) }
  d <- matched[!is.na(stage_group), c("stage_group", sc), with = FALSE]
  dl <- melt(d, id.vars = "stage_group", measure.vars = sc,
             variable.name = "sig", value.name = "score")
  dl <- dl[is.finite(score)]
  dl[, sig := gsub("_", "-", sub("^score_", "", sig))]   # DE_only -> DE-only
  pl <- copy(stage_tab); setnames(pl, "signature", "sig", skip_absent = TRUE)
  pl <- pl[sig %in% unique(dl$sig)]
  pl[, plab := sprintf("ANOVA p = %.3f", anova_p)]
  STG <- c("Stage I" = "#4DBBD5", "Stage II" = "#7CC68E",
           "Stage III" = "#E07B5E", "Stage IV" = "#B40426")
  p <- ggplot(dl, aes(x = stage_group, y = score, fill = stage_group)) +
    ggdist::stat_halfeye(adjust = 0.8, width = 0.55, .width = 0,
                         justification = -0.12, point_colour = NA, alpha = 0.55) +
    geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 0.9,
                 color = "grey25", linewidth = 0.3) +
    ggbeeswarm::geom_quasirandom(width = 0.09, size = 0.25, alpha = 0.16,
                                 shape = 16) +
    geom_text(data = pl, aes(x = 0.7, y = Inf, label = plab),
              inherit.aes = FALSE, hjust = 0, vjust = 1.6, size = 2.4,
              color = "grey25", fontface = "italic") +
    facet_wrap(~sig, scales = "free_y") +
    scale_fill_manual(values = STG, name = "AJCC stage") +
    labs(title = "**SGCB channel signatures by AJCC stage (TCGA-KIRC)**",
         subtitle = "Per-patient signature score across pathologic stage; complements the Cox/KM panels of Fig 7 (no overlap)",
         x = NULL, y = "Signature score (z)") +
    THEME_JBHI(base_size = 10) +
    theme(plot.title = ggtext::element_markdown(size = 11, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 8, color = "grey40"),
          strip.text = element_text(face = "bold", size = 9),
          axis.text.x = element_text(size = 8),
          legend.position = "bottom")
  SAVE_FIG_V5("FigS5_tcga_kirc_clinical", p, 9.5, 5.0)
}

# --- run all redesigns (guarded) -----------------------------------------
.try_redesign <- function(nm, fn) {
  tryCatch({ fn(); cat("  [37d] redesigned", nm, "\n") },
           error = function(e) cat("  [37d][ERROR]", nm, ":",
                                   conditionMessage(e), "\n"))
}
.try_redesign("Fig2", redesign_fig2)
.try_redesign("Fig3", redesign_fig3)
.try_redesign("Fig4", redesign_fig4)
.try_redesign("Fig6", redesign_fig6)
.try_redesign("Fig7_seqc", redesign_fig7_seqc)
.try_redesign("Fig7", redesign_fig7)
.try_redesign("Fig8", redesign_fig8)
.try_redesign("FigS1", redesign_figS1)
.try_redesign("FigS3", redesign_figS3)
.try_redesign("FigS5", redesign_figS5)
cat("=== 37d overlay done ===\n")
