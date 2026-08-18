# =============================================================================
#
# Ablates the 5 binary switches exposed by sgcbDE():
#
# For each (switch, dataset, seed) we run sgcbDE() once, and record DE
# detection quality at FDR 5% against the simulation truth plus wall-clock
# time. Two sample-size regimes are included: n=5/group triggers the
# small-sample regularizers (Firth, gamma_submodel); n=10/group probes the
# regime where small-sample guards are less active so we can separate their
# effect from the other switches.
#
# Output: output/ablation_sgcb_metrics.csv
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
    library(SGCB)
    library(data.table)
})

DATA_DIR    <- paste0(PROJECT_ROOT, "/benchmark/data")
RESULTS_DIR <- paste0(PROJECT_ROOT, "/benchmark/output")
dir.create(DATA_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# NB simulation: same generator used in the main benchmark (19_bulk_benchmark.R
# + _archive/generate_sim_data.R). Kept local so the ablation script is
# self-contained and can be re-run without external state.
# -----------------------------------------------------------------------------
generate_sim_data <- function(n_genes = 5000L, n_per_group = 5L,
                              pDE = 0.1, effect_size = 2.0, seed = 1L) {
    set.seed(seed)
    n_samples <- n_per_group * 2L
    base_mean  <- pmax(exp(rnorm(n_genes, mean = 4, sd = 2)), 1)
    dispersion <- 0.1 + 1 / sqrt(base_mean)
    n_de       <- round(n_genes * pDE)
    de_idx     <- sample.int(n_genes, n_de)
    de_dir     <- sample(c(-1, 1), n_de, replace = TRUE)
    counts <- matrix(0L, n_genes, n_samples,
                     dimnames = list(paste0("gene_", seq_len(n_genes)),
                                     paste0("sample_", seq_len(n_samples))))
    for (i in seq_len(n_genes)) {
        mu   <- base_mean[i]
        size <- 1 / dispersion[i]
        counts[i, 1:n_per_group] <- rnbinom(n_per_group, size = size, mu = mu)
        mu_t <- mu
        if (i %in% de_idx) {
            d    <- de_dir[match(i, de_idx)]
            mu_t <- if (d == 1L) mu * effect_size else mu / effect_size
        }
        counts[i, (n_per_group + 1):n_samples] <- rnbinom(n_per_group, size = size, mu = mu_t)
    }
    is_de <- logical(n_genes); is_de[de_idx] <- TRUE
    list(counts = counts,
         truth  = data.table(gene = rownames(counts), is_de = is_de),
         group  = factor(c(rep("ctrl", n_per_group), rep("treat", n_per_group))))
}

# -----------------------------------------------------------------------------
# Ablation grid. `base` has all switches on. Each named variant flips exactly
# one switch to FALSE so every row isolates one component's contribution.
# -----------------------------------------------------------------------------
variant_grid <- list(
    base                = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_natural_grad     = list(use_natural_grad = FALSE, use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_hierarchical     = list(use_natural_grad = TRUE,  use_hierarchical = FALSE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_firth            = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = FALSE, use_gamma_submodel = TRUE,
                               use_gg_variance  = TRUE),
    no_gamma_submodel   = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = FALSE,
                               use_gg_variance  = TRUE),
    no_gg_variance      = list(use_natural_grad = TRUE,  use_hierarchical = TRUE,
                               use_firth        = TRUE,  use_gamma_submodel = TRUE,
                               use_gg_variance  = FALSE)
)

# -----------------------------------------------------------------------------
# Per-run wrapper. sgcbDE() filters low-count genes (min_count = 10) so its
# row order differs from sim$counts; we therefore align by gene name rather
# than by position (position-based comparison produced FDR > 1 because R
# recycled length-mismatched logical vectors).
#
# The NB simulation injects mean shifts only, so:
#   - DE ground truth is the sim$truth$is_de vector (for F1 / actual FDR);
#   - every DV and DG detection at nominal 5% is by construction a null false
#     positive, which lets us read off DV/DG null FPR per ablation variant.
# -----------------------------------------------------------------------------
run_one <- function(sim, switches, alpha = 0.05) {
    args <- c(list(counts = sim$counts, group = sim$group), switches)
    t0   <- Sys.time()
    res  <- tryCatch(do.call(sgcbDE, args), error = function(e) e)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
        return(data.table(F1 = NA_real_, de_FDR = NA_real_, de_n_sig = NA_integer_,
                          dv_FPR = NA_real_, dv_n_sig = NA_integer_,
                          dg_FPR = NA_real_, dg_n_sig = NA_integer_,
                          n_tested = NA_integer_, time_sec = tsec,
                          error = conditionMessage(res)))
    }

    genes_tested <- rownames(res)
    truth_lookup <- setNames(sim$truth$is_de, sim$truth$gene)
    truth_tested <- truth_lookup[genes_tested]

    # -- DE channel vs ground truth -------------------------------------------
    padj_de <- res$padj
    sig_de  <- !is.na(padj_de) & padj_de < alpha
    TP <- sum(sig_de &  truth_tested)
    FP <- sum(sig_de & !truth_tested)
    FN <- sum(!sig_de &  truth_tested)
    prec <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
    rec  <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
    F1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
                2 * prec * rec / (prec + rec) else NA_real_
    de_FDR <- if (sum(sig_de) > 0) FP / sum(sig_de) else 0

    # -- DV channel null FPR (sim has no DV signal) ---------------------------
    padj_dv <- res$dv_padj
    sig_dv  <- !is.na(padj_dv) & padj_dv < alpha
    dv_FPR  <- mean(sig_dv)

    # -- DG channel null FPR (sim has no DG signal) ---------------------------
    padj_dg <- res$dg_shape_padj
    sig_dg  <- !is.na(padj_dg) & padj_dg < alpha
    dg_FPR  <- mean(sig_dg)

    data.table(F1         = F1,
               de_FDR     = de_FDR,
               de_n_sig   = sum(sig_de),
               dv_FPR     = dv_FPR,
               dv_n_sig   = sum(sig_dv),
               dg_FPR     = dg_FPR,
               dg_n_sig   = sum(sig_dg),
               n_tested   = length(genes_tested),
               time_sec   = tsec,
               error      = NA_character_)
}

# -----------------------------------------------------------------------------
# Drive the grid. We use 3 seeds per (n_per_group, effect_size) and average.
# -----------------------------------------------------------------------------
cfgs  <- expand.grid(n_per_group = c(5L, 10L),
                     effect_size = 2.0,
                     seed        = 1:3,
                     KEEP.OUT.ATTRS = FALSE)
rows  <- list()
total <- nrow(cfgs) * length(variant_grid)
k <- 0L

for (r in seq_len(nrow(cfgs))) {
    sim <- generate_sim_data(n_genes     = 5000L,
                             n_per_group = cfgs$n_per_group[r],
                             pDE         = 0.1,
                             effect_size = cfgs$effect_size[r],
                             seed        = cfgs$seed[r])
    for (vname in names(variant_grid)) {
        k <- k + 1L
        cat(sprintf("[%2d/%2d] n=%d ef=%.1f seed=%d | %-18s ... ",
                    k, total, cfgs$n_per_group[r], cfgs$effect_size[r],
                    cfgs$seed[r], vname))
        out <- run_one(sim, variant_grid[[vname]])
        cat(sprintf("F1=%.3f deFDR=%.3f dvFPR=%.4f dgFPR=%.4f t=%.1fs%s\n",
                    out$F1, out$de_FDR, out$dv_FPR, out$dg_FPR, out$time_sec,
                    if (!is.na(out$error)) paste0(" ERR:", out$error) else ""))
        rows[[length(rows) + 1L]] <- cbind(
            data.table(variant     = vname,
                       n_per_group = cfgs$n_per_group[r],
                       effect_size = cfgs$effect_size[r],
                       seed        = cfgs$seed[r]),
            out)
    }
}

all_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
out_csv <- file.path(RESULTS_DIR, "ablation_sgcb_metrics.csv")
fwrite(all_dt, out_csv)
cat("\nSaved raw metrics:", out_csv, "\n")

# -----------------------------------------------------------------------------
# Aggregate summary: mean over seeds per (variant, n). Delta-vs-base added for
# quick reading in the console summary.
# -----------------------------------------------------------------------------
agg <- all_dt[!is.na(F1), .(mean_F1    = mean(F1),
                             mean_deFDR = mean(de_FDR),
                             mean_dvFPR = mean(dv_FPR),
                             mean_dgFPR = mean(dg_FPR),
                             mean_deSig = mean(de_n_sig),
                             mean_dvSig = mean(dv_n_sig),
                             mean_dgSig = mean(dg_n_sig),
                             mean_time  = mean(time_sec),
                             sd_F1      = sd(F1),
                             sd_deFDR   = sd(de_FDR),
                             sd_dvFPR   = sd(dv_FPR)),
                  by = .(variant, n_per_group)]

for (n_now in unique(agg$n_per_group)) {
    b <- agg[variant == "base" & n_per_group == n_now]
    agg[n_per_group == n_now,
        `:=`(delta_F1    = mean_F1    - b$mean_F1,
             delta_deFDR = mean_deFDR - b$mean_deFDR,
             delta_dvFPR = mean_dvFPR - b$mean_dvFPR,
             delta_dgFPR = mean_dgFPR - b$mean_dgFPR,
             time_ratio  = mean_time  / b$mean_time)]
}

out_csv2 <- file.path(RESULTS_DIR, "ablation_sgcb_summary.csv")
fwrite(agg, out_csv2)
cat("Saved summary:     ", out_csv2, "\n\n")

cat("=== Summary ===\n")
print(agg[order(n_per_group, variant)])
cat("\nDone.\n")
