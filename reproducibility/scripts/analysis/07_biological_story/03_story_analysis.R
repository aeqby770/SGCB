PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
suppressPackageStartupMessages({
  library(data.table)
})

source(paste0(PROJECT_ROOT, "/biological_story_workspace/00_story_config.R"))
source(paste0(PROJECT_ROOT, "/biological_story_workspace/01_story_utils.R"))

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("STORY_DEFAULT_CASES", "signature_score", "channel_category", "variable", "variable_level"))
}

story_parse_case_ids <- function(args) {
  case_arg <- if (length(args) >= 1L) args[1] else NA_character_
  if (is.na(case_arg) || identical(case_arg, "") || identical(case_arg, "all")) {
    return(STORY_DEFAULT_CASES)
  }
  strsplit(case_arg, ",", fixed = TRUE)[[1]]
}

story_parse_flag <- function(args, flag) {
  any(args[-1] %in% flag)
}

story_safe_test <- function(x, g, method = c("wilcox", "kruskal")) {
  method <- match.arg(method)
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  if (length(x) == 0L || length(unique(g)) < 2L) {
    return(NA_real_)
  }
  out <- tryCatch(
    if (method == "wilcox") {
      stats::wilcox.test(x ~ g)$p.value
    } else {
      stats::kruskal.test(x ~ g)$p.value
    },
    error = function(e) NA_real_
  )
  as.numeric(out)
}

story_clinical_associations <- function(sample_scores) {
  if (nrow(sample_scores) == 0L) {
    return(data.table())
  }
  vars <- intersect(c("group", "tumor_grade", "tumor_stage", "story_group", "primary_diagnosis"), names(sample_scores))
  vars <- vars[vars != "channel_category"]
  if (length(vars) == 0L) {
    return(data.table())
  }
  assoc_list <- lapply(vars, function(vn) {
    dt <- as.data.table(sample_scores[!is.na(sample_scores[[vn]]) & !is.na(sample_scores[["signature_score"]]), , drop = FALSE])
    if (nrow(dt) == 0L) {
      return(data.table())
    }
    split_levels <- unique(dt[[vn]])
    stat_method <- if (length(split_levels) == 2L) "wilcox" else "kruskal"
    channel_vals <- unique(dt$channel_category)
    out_rows <- lapply(channel_vals, function(ch) {
      sub_dt <- dt[channel_category == ch]
      if (nrow(sub_dt) == 0L) {
        return(data.table())
      }
      pval <- story_safe_test(sub_dt[["signature_score"]], sub_dt[[vn]], method = stat_method)
      med_dt <- as.data.table(sub_dt)[, list(
        median_signature_score = median(signature_score, na.rm = TRUE),
        mean_signature_score = mean(signature_score, na.rm = TRUE),
        n = .N
      ), by = vn]
      setnames(med_dt, old = vn, new = "variable_level")
      med_dt[, channel_category := ch]
      med_dt[, variable := vn]
      med_dt[, p_value := pval]
      setcolorder(med_dt, c("channel_category", "variable", "variable_level", "p_value", "median_signature_score", "mean_signature_score", "n"))
      med_dt
    })
    rbindlist(out_rows, use.names = TRUE, fill = TRUE)
  })
  rbindlist(assoc_list, use.names = TRUE, fill = TRUE)
}

story_run_case_analysis <- function(case_id, rerun = FALSE, use_cached = TRUE) {
  case_cfg <- story_get_case(case_id)
  dirs <- story_case_dirs(case_cfg)
  case_data <- story_read_case_data(case_cfg)
  result_dt <- story_load_or_run_result(case_cfg, case_data, use_cached = use_cached, rerun = rerun)
  standardized_dt <- story_standardize_result(result_dt, fdr_cutoff = STORY_FDR_CUTOFF)
  channel_summary_dt <- story_channel_summary(standardized_dt)
  representative_dt <- story_pick_representative_genes(standardized_dt, top_n = STORY_TOP_FEATURES_PER_CHANNEL)
  enrichment_dt <- story_run_case_enrichment(standardized_dt)
  marker_ref_dt <- story_case_marker_reference(case_cfg, case_data)
  marker_enrichment_dt <- story_marker_enrichment(standardized_dt, marker_ref_dt)
  representative_marker_dt <- story_representative_marker_map(representative_dt, marker_ref_dt)
  sample_score_dt <- story_compute_sample_scores(
    log_expr = case_data$log_expr,
    standardized_dt = standardized_dt,
    sample_ids = case_data$sample_ids,
    group = case_data$group,
    clinical = case_data$clinical,
    top_k = STORY_SIGNATURE_TOP_K
  )
  marker_score_dt <- story_marker_scores(
    log_expr = case_data$log_expr,
    sample_ids = case_data$sample_ids,
    group = case_data$group,
    marker_ref_dt = marker_ref_dt,
    clinical = case_data$clinical
  )
  distribution_dt <- story_prepare_distribution_data(
    log_expr = case_data$log_expr,
    sample_ids = case_data$sample_ids,
    group = case_data$group,
    representative_dt = representative_dt,
    clinical = case_data$clinical
  )
  clinical_assoc_dt <- story_clinical_associations(sample_score_dt)
  sgcb_result_source <- if (isTRUE(use_cached && !rerun && !is.null(case_cfg$cached_result_path))) {
    case_cfg$cached_result_path
  } else {
    file.path(dirs$cache, paste0(case_cfg$case_id, "_sgcb_result.rds"))
  }
  case_meta_dt <- data.table(
    case_id = case_cfg$case_id,
    title = case_cfg$title,
    modality = case_cfg$modality,
    n_features = nrow(case_data$expr),
    n_samples = ncol(case_data$expr),
    group_levels = paste(levels(case_data$group), collapse = ";"),
    group_counts = paste(as.integer(table(case_data$group)), collapse = ";"),
    cell_types = if (!is.null(case_data$cell_types)) paste(stats::na.omit(case_data$cell_types), collapse = ";") else NA_character_,
    used_cached_result = isTRUE(use_cached && !rerun && !is.null(case_cfg$cached_result_path)),
    sgcb_result_source = sgcb_result_source
  )
  saveRDS(result_dt, file.path(dirs$cache, paste0(case_cfg$case_id, "_sgcb_result.rds")))
  saveRDS(case_data, file.path(dirs$cache, paste0(case_cfg$case_id, "_case_data.rds")))
  fwrite(case_meta_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_case_metadata.csv")))
  fwrite(standardized_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_sgcb_channels.csv")))
  fwrite(channel_summary_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_channel_summary.csv")))
  fwrite(representative_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_representative_features.csv")))
  fwrite(enrichment_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_enrichment.csv")))
  fwrite(marker_ref_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_celltype_marker_reference.csv")))
  fwrite(marker_enrichment_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_celltype_marker_enrichment.csv")))
  fwrite(representative_marker_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_representative_marker_map.csv")))
  fwrite(sample_score_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_sample_signature_scores.csv")))
  fwrite(marker_score_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_celltype_marker_scores.csv")))
  fwrite(distribution_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_representative_distributions_long.csv")))
  fwrite(clinical_assoc_dt, file.path(dirs$tables, paste0(case_cfg$case_id, "_clinical_associations.csv")))
  invisible(list(
    case_cfg = case_cfg,
    dirs = dirs,
    case_data = case_data,
    result_dt = result_dt,
    standardized_dt = standardized_dt,
    channel_summary_dt = channel_summary_dt,
    representative_dt = representative_dt,
    enrichment_dt = enrichment_dt,
    marker_ref_dt = marker_ref_dt,
    marker_enrichment_dt = marker_enrichment_dt,
    representative_marker_dt = representative_marker_dt,
    sample_score_dt = sample_score_dt,
    marker_score_dt = marker_score_dt,
    distribution_dt = distribution_dt,
    clinical_assoc_dt = clinical_assoc_dt
  ))
}

story_run_many_cases <- function(case_ids, rerun = FALSE, use_cached = TRUE) {
  outputs <- lapply(case_ids, function(case_id) {
    message("[story] analyzing ", case_id)
    story_run_case_analysis(case_id = case_id, rerun = rerun, use_cached = use_cached)
  })
  names(outputs) <- case_ids
  invisible(outputs)
}

story_analysis_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  case_ids <- story_parse_case_ids(args)
  rerun <- story_parse_flag(args, c("--rerun", "rerun"))
  no_cache <- story_parse_flag(args, c("--no-cache", "no-cache"))
  story_run_many_cases(case_ids = case_ids, rerun = rerun, use_cached = !no_cache)
}

if (sys.nframe() == 0L) {
  story_analysis_main()
}
