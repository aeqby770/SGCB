# =============================================================================
#
# Complements 50_ablation_sgcb.R (mean-shift only) with two additional regimes
# that actually exercise the GG-side components:
#
#                     (probes Firth / hierarchical prior / gamma submodel under
#                      sparse-data stress; we track p-value NA rate per gene as a
#                      proxy for convergence failures.)
#                     so the DG channel has a real signal to recover.
#
# For each (scenario, variant, seed) run we record DE F1/FDR, DV null FPR, DG
# FPR, DG recall (scenario C), p-value NA rate, and wall-clock time. Results
# are written to ablation_sgcb_extended_{metrics,summary}.csv.
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
# Scenario A: classic mean-shift NB (same as 50_ablation_sgcb.R, kept for cross-
# scenario comparability with the new regimes).
# -----------------------------------------------------------------------------
gen_mean_shift <- function(n_genes = 5000L, n_per_group = 5L, pDE = 0.1,
                           effect_size = 2.0, seed = 1L) {
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
         truth  = data.table(gene = rownames(counts), is_de = is_de,
                             is_dg = FALSE),
         group  = factor(c(rep("ctrl", n_per_group),
                           rep("treat", n_per_group))))
}

# -----------------------------------------------------------------------------
# Scenario B: small-sample stress. n=3/group, low-count baseline (rnorm mean=2
# instead of 4) and inflated dispersion (0.2 + 1.5/sqrt(mu) instead of the
# benchmark default). These are the regimes where the Firth penalty and the
# hierarchical prior are motivated in the Methods.
# -----------------------------------------------------------------------------
gen_stress <- function(n_genes = 3000L, n_per_group = 3L, pDE = 0.1,
                       effect_size = 2.0, seed = 1L) {
    set.seed(seed)
    n_samples <- n_per_group * 2L
    base_mean  <- pmax(exp(rnorm(n_genes, mean = 2, sd = 1.3)), 1)
    dispersion <- 0.2 + 1.5 / sqrt(base_mean)
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
         truth  = data.table(gene = rownames(counts), is_de = is_de,
                             is_dg = FALSE),
         group  = factor(c(rep("ctrl", n_per_group),
                           rep("treat", n_per_group))))
}

# -----------------------------------------------------------------------------
# Scenario C: shape-shift simulation. For a fraction pDG of genes we replace
# group-2 counts with a mean-matched bimodal mixture: half of the samples are
# drawn from NB(0.2*mu, phi), the other half from NB(1.8*mu, phi). The two
# components give E[X] = mu exactly, but the resulting per-gene distribution
# DE/DV signals are NOT injected, so:
#   DE calls on non-DG genes  = null  (FDR target: <= alpha)
#   DG calls on DG-truth genes = real  (recall target: > 0)
#   DG calls on non-DG genes  = null  (FDR on DG channel)
# -----------------------------------------------------------------------------
gen_shape_shift <- function(n_genes = 3000L, n_per_group = 10L, pDG = 0.2,
                            seed = 1L) {
    set.seed(seed)
    if (n_per_group %% 2L != 0L)
        stop("n_per_group must be even for the 50/50 mixture.")
    n_samples <- n_per_group * 2L
    base_mean  <- pmax(exp(rnorm(n_genes, mean = 4, sd = 2)), 1)
    dispersion <- 0.1 + 1 / sqrt(base_mean)
    n_dg       <- round(n_genes * pDG)
    dg_idx     <- sample.int(n_genes, n_dg)
    counts <- matrix(0L, n_genes, n_samples,
                     dimnames = list(paste0("gene_", seq_len(n_genes)),
                                     paste0("sample_", seq_len(n_samples))))
    half <- n_per_group %/% 2L
    for (i in seq_len(n_genes)) {
        mu   <- base_mean[i]
        size <- 1 / dispersion[i]
        # control arm: straight NB
        counts[i, 1:n_per_group] <- rnbinom(n_per_group, size = size, mu = mu)
        # treatment arm: NB for non-DG, bimodal mixture for DG
        if (i %in% dg_idx) {
            lo <- rnbinom(half, size = size, mu = 0.2 * mu)
            hi <- rnbinom(half, size = size, mu = 1.8 * mu)
            counts[i, (n_per_group + 1):n_samples] <- sample(c(lo, hi))
        } else {
            counts[i, (n_per_group + 1):n_samples] <-
                rnbinom(n_per_group, size = size, mu = mu)
        }
    }
    is_dg <- logical(n_genes); is_dg[dg_idx] <- TRUE
    list(counts = counts,
         truth  = data.table(gene = rownames(counts),
                             is_de = FALSE,
                             is_dg = is_dg),
         group  = factor(c(rep("ctrl", n_per_group),
                           rep("treat", n_per_group))))
}

# -----------------------------------------------------------------------------
# Ablation grid: one switch off at a time.
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
# Run wrapper. Align by gene name (sgcbDE filters low-count genes so row order
# differs from sim$counts). Also reports:
#     values flag optimization/convergence trouble on sparse data.
# -----------------------------------------------------------------------------
run_one <- function(sim, switches, alpha = 0.05) {
    args <- c(list(counts = sim$counts, group = sim$group), switches)
    t0   <- Sys.time()
    res  <- tryCatch(do.call(sgcbDE, args), error = function(e) e)
    tsec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
        return(data.table(F1 = NA_real_, de_FDR = NA_real_, de_n_sig = NA_integer_,
                          de_na_rate = 1,
                          dv_FPR = NA_real_, dv_n_sig = NA_integer_,
                          dg_FPR = NA_real_, dg_n_sig = NA_integer_,
                          dg_recall = NA_real_,
                          n_tested = 0L, time_sec = tsec,
                          error = conditionMessage(res)))
    }

    genes_tested <- rownames(res)
    truth_de <- setNames(sim$truth$is_de, sim$truth$gene)[genes_tested]
    truth_dg <- setNames(sim$truth$is_dg, sim$truth$gene)[genes_tested]

    # DE
    padj_de    <- res$padj
    sig_de     <- !is.na(padj_de) & padj_de < alpha
    de_na_rate <- mean(is.na(padj_de))
    TP <- sum(sig_de &  truth_de)
    FP <- sum(sig_de & !truth_de)
    FN <- sum(!sig_de &  truth_de)
    prec <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
    rec  <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
    F1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
                2 * prec * rec / (prec + rec) else NA_real_
    de_FDR <- if (sum(sig_de) > 0) FP / sum(sig_de) else 0

    # null-FPR on complement. When truth set is empty (scenarios A/B), the
    # "null FPR" is just the overall call rate and recall is NA.
    channel_stats <- function(sig, truth) {
        if (any(truth)) {
            rec  <- sum(sig &  truth) / sum(truth)
            nullFPR <- sum(sig & !truth) / max(sum(!truth), 1L)
        } else {
            rec <- NA_real_
            nullFPR <- mean(sig)
        }
        list(recall = rec, null_FPR = nullFPR, n_sig = sum(sig))
    }

    padj_dv <- res$dv_padj
    sig_dv  <- !is.na(padj_dv) & padj_dv < alpha
    dv <- channel_stats(sig_dv, truth_dg)

    padj_dg <- res$dg_shape_padj
    sig_dg  <- !is.na(padj_dg) & padj_dg < alpha
    dg <- channel_stats(sig_dg, truth_dg)

    data.table(F1         = F1,
               de_FDR     = de_FDR,
               de_n_sig   = sum(sig_de),
               de_na_rate = de_na_rate,
               dv_recall  = dv$recall,
               dv_FPR     = dv$null_FPR,
               dv_n_sig   = dv$n_sig,
               dg_recall  = dg$recall,
               dg_FPR     = dg$null_FPR,
               dg_n_sig   = dg$n_sig,
               n_tested   = length(genes_tested),
               time_sec   = tsec,
               error      = NA_character_)
}

# -----------------------------------------------------------------------------
# Scenario drive table. Each row spawns one sim, then we run every variant.
# -----------------------------------------------------------------------------
scenarios <- rbindlist(list(
    data.table(scenario = "mean_shift", n_per_group = 5L,  seed = 1:3),
    data.table(scenario = "mean_shift", n_per_group = 10L, seed = 1:3),
    data.table(scenario = "stress",     n_per_group = 3L,  seed = 1:3),
    data.table(scenario = "shape_shift",n_per_group = 10L, seed = 1:3)
))

rows  <- list()
total <- nrow(scenarios) * length(variant_grid)
k <- 0L

for (r in seq_len(nrow(scenarios))) {
    sc <- scenarios$scenario[r]
    n  <- scenarios$n_per_group[r]
    sd <- scenarios$seed[r]
    sim <- switch(sc,
        mean_shift  = gen_mean_shift (n_per_group = n, seed = sd),
        stress      = gen_stress     (n_per_group = n, seed = sd),
        shape_shift = gen_shape_shift(n_per_group = n, seed = sd)
    )
    for (vname in names(variant_grid)) {
        k <- k + 1L
        cat(sprintf("[%2d/%2d] %-12s n=%2d seed=%d | %-18s ... ",
                    k, total, sc, n, sd, vname))
        out <- run_one(sim, variant_grid[[vname]])
        cat(sprintf("F1=%.3f deFDR=%.3f dvFPR=%.4f dvRec=%s dgRec=%s t=%.1fs%s\n",
                    out$F1, out$de_FDR, out$dv_FPR,
                    if (is.na(out$dv_recall)) "  NA" else sprintf("%.3f", out$dv_recall),
                    if (is.na(out$dg_recall)) "  NA" else sprintf("%.3f", out$dg_recall),
                    out$time_sec,
                    if (!is.na(out$error)) paste0(" ERR:", out$error) else ""))
        rows[[length(rows) + 1L]] <- cbind(
            data.table(scenario    = sc,
                       variant     = vname,
                       n_per_group = n,
                       seed        = sd),
            out)
    }
}

all_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
out_csv <- file.path(RESULTS_DIR, "ablation_sgcb_extended_metrics.csv")
fwrite(all_dt, out_csv)
cat("\nSaved raw metrics:", out_csv, "\n")

# -----------------------------------------------------------------------------
# Per (scenario, variant) means, plus delta-vs-base and time ratio.
# -----------------------------------------------------------------------------
agg <- all_dt[, .(mean_F1      = mean(F1,         na.rm = TRUE),
                  mean_deFDR   = mean(de_FDR,     na.rm = TRUE),
                  mean_deNA    = mean(de_na_rate, na.rm = TRUE),
                  mean_dvRec   = mean(dv_recall,  na.rm = TRUE),
                  mean_dvFPR   = mean(dv_FPR,     na.rm = TRUE),
                  mean_dgRec   = mean(dg_recall,  na.rm = TRUE),
                  mean_dgFPR   = mean(dg_FPR,     na.rm = TRUE),
                  mean_deSig   = mean(de_n_sig,   na.rm = TRUE),
                  mean_dvSig   = mean(dv_n_sig,   na.rm = TRUE),
                  mean_dgSig   = mean(dg_n_sig,   na.rm = TRUE),
                  mean_time    = mean(time_sec,   na.rm = TRUE),
                  n_ok         = sum(!is.na(F1))),
              by = .(scenario, variant, n_per_group)]

for (sc in unique(agg$scenario)) {
    for (n_now in unique(agg[scenario == sc]$n_per_group)) {
        b <- agg[scenario == sc & variant == "base" & n_per_group == n_now]
        if (nrow(b) == 0L) next
        agg[scenario == sc & n_per_group == n_now,
            `:=`(delta_F1     = mean_F1     - b$mean_F1,
                 delta_deFDR  = mean_deFDR  - b$mean_deFDR,
                 delta_dvRec  = mean_dvRec  - b$mean_dvRec,
                 delta_dvFPR  = mean_dvFPR  - b$mean_dvFPR,
                 delta_dgRec  = mean_dgRec  - b$mean_dgRec,
                 delta_dgFPR  = mean_dgFPR  - b$mean_dgFPR,
                 delta_deNA   = mean_deNA   - b$mean_deNA,
                 time_ratio   = mean_time   / b$mean_time)]
    }
}

out_csv2 <- file.path(RESULTS_DIR, "ablation_sgcb_extended_summary.csv")
fwrite(agg, out_csv2)
cat("Saved summary:     ", out_csv2, "\n\n")

cat("=== Summary (per scenario) ===\n")
for (sc in unique(agg$scenario)) {
    cat("\n-- scenario:", sc, "--\n")
    print(agg[scenario == sc][order(n_per_group, variant)])
}
cat("\nDone.\n")
