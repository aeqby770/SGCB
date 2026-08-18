suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(SGCB)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

story_annotation_cache_dir <- function() {
  root_dir <- get0("STORY_ROOT", ifnotfound = normalizePath(".", winslash = "/", mustWork = TRUE))
  cache_dir <- file.path(root_dir, "cache", "annotations")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

story_local_go_maps <- function() {
  cache_file <- file.path(story_annotation_cache_dir(), "go_term_maps_human.rds")
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  entrez_ids <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")
  go_map <- as.data.table(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_ids,
    columns = c("GO", "ONTOLOGY"),
    keytype = "ENTREZID"
  ))
  go_map <- unique(go_map[!is.na(GO) & !is.na(ONTOLOGY), .(ENTREZID, GO, ONTOLOGY)])
  go_name <- as.data.table(AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(go_map$GO),
    columns = "TERM",
    keytype = "GOID"
  ))
  go_name <- unique(go_name[!is.na(TERM), .(GOID, TERM)])
  out <- list(
    term2gene = go_map,
    term2name = go_name
  )
  saveRDS(out, cache_file)
  out
}

story_download_text_with_curl <- function(url, dest_file) {
  curl_path <- Sys.which("curl.exe")
  if (!nzchar(curl_path)) {
    return(FALSE)
  }
  txt <- tryCatch(
    suppressWarnings(system2(
      curl_path,
      args = c("-k", "-L", "-sS", "--connect-timeout", "60", "--max-time", "180", url),
      stdout = TRUE,
      stderr = FALSE
    )),
    error = function(e) character()
  )
  txt <- txt[nzchar(txt)]
  if (length(txt) == 0L) {
    return(FALSE)
  }
  writeLines(enc2utf8(txt), dest_file, useBytes = TRUE)
  file.exists(dest_file) && file.info(dest_file)$size > 0L
}

story_download_kegg_tables <- function() {
  link_url <- "http://rest.kegg.jp/link/pathway/hsa"
  list_url <- "http://rest.kegg.jp/list/pathway/hsa"
  link_file <- file.path(story_annotation_cache_dir(), "kegg_link_pathway_hsa.tsv")
  list_file <- file.path(story_annotation_cache_dir(), "kegg_list_pathway_hsa.tsv")
  if (!file.exists(link_file) || file.info(link_file)$size <= 0L) {
    ok <- story_download_text_with_curl(link_url, link_file)
    if (!ok) {
      link_txt <- readLines(link_url, warn = FALSE, encoding = "UTF-8")
      writeLines(link_txt, link_file, useBytes = TRUE)
    }
  }
  if (!file.exists(list_file) || file.info(list_file)$size <= 0L) {
    ok <- story_download_text_with_curl(list_url, list_file)
    if (!ok) {
      list_txt <- readLines(list_url, warn = FALSE, encoding = "UTF-8")
      writeLines(list_txt, list_file, useBytes = TRUE)
    }
  }
  link_txt <- readLines(link_file, warn = FALSE, encoding = "UTF-8")
  list_txt <- readLines(list_file, warn = FALSE, encoding = "UTF-8")
  list(link = link_txt, list = list_txt)
}

story_local_kegg_maps <- function(refresh = FALSE) {
  cache_file <- file.path(story_annotation_cache_dir(), "kegg_term_maps_hsa.rds")
  if (!refresh && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  kegg_raw <- story_download_kegg_tables()
  link_dt <- fread(
    text = paste(kegg_raw$link, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    col.names = c("gene_id", "pathway_id"),
    encoding = "UTF-8"
  )
  link_dt[, ENTREZID := sub("^hsa:", "", gene_id)]
  link_dt[, pathway_id := sub("^path:", "", pathway_id)]
  link_dt <- unique(link_dt[, .(pathway_id, ENTREZID)])
  name_dt <- fread(
    text = paste(kegg_raw$list, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    col.names = c("pathway_id", "pathway_name"),
    encoding = "UTF-8"
  )
  name_dt[, pathway_name := sub(" - Homo sapiens [(]human[)]$", "", pathway_name)]
  name_dt <- unique(name_dt[, .(pathway_id, pathway_name)])
  out <- list(
    term2gene = link_dt,
    term2name = name_dt
  )
  saveRDS(out, cache_file)
  out
}

story_local_enricher <- function(gene_entrez, universe_entrez, term2gene, term2name) {
  gene_entrez <- unique(as.character(gene_entrez))
  universe_entrez <- unique(as.character(universe_entrez))
  term2gene <- as.data.table(copy(term2gene))
  term2name <- as.data.table(copy(term2name))
  term_col <- names(term2gene)[1]
  gene_col <- names(term2gene)[2]
  name_col <- names(term2name)[2]
  term2gene <- unique(term2gene[get(gene_col) %in% universe_entrez])
  if (nrow(term2gene) == 0L) {
    return(NULL)
  }
  term2name <- unique(term2name[get(names(term2name)[1]) %in% unique(term2gene[[term_col]])])
  suppressMessages(clusterProfiler::enricher(
    gene = gene_entrez,
    universe = universe_entrez,
    TERM2GENE = term2gene[, .(term = get(term_col), gene = get(gene_col))],
    TERM2NAME = term2name[, .(term = get(names(term2name)[1]), name = get(name_col))],
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  ))
}

story_case_dirs <- function(case_cfg) {
  dirs <- list(
    root = case_cfg$output_dir,
    tables = file.path(case_cfg$output_dir, "tables"),
    figures = file.path(case_cfg$output_dir, "figures"),
    cache = file.path(case_cfg$output_dir, "cache")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

story_read_bulk_case <- function(case_cfg) {
  dat <- readRDS(case_cfg$input_path)
  counts <- as.matrix(dat$counts)
  rownames(counts) <- make.unique(rownames(counts))
  group <- factor(dat$group)
  if (!is.null(case_cfg$group_levels)) {
    group <- factor(group, levels = case_cfg$group_levels)
  }
  clinical <- NULL
  if ("metadata" %in% names(dat) && !is.null(dat$metadata)) {
    clinical <- as.data.table(copy(dat$metadata))
    if (!"sample_id" %in% names(clinical)) {
      if ("sample" %in% names(clinical)) {
        setnames(clinical, old = "sample", new = "sample_id")
      } else {
        clinical[, sample_id := colnames(counts)]
      }
    }
    clinical <- clinical[match(colnames(counts), sample_id)]
  }
  list(
    expr = counts,
    log_expr = edgeR::cpm(counts, log = TRUE, prior.count = 1),
    group = droplevels(group),
    sample_ids = colnames(counts),
    clinical = clinical,
    cell_types = if ("cell_types" %in% names(dat)) as.character(dat$cell_types) else NULL,
    feature_key = "gene",
    assay_name = "counts"
  )
}

story_read_proteomics_case <- function(case_cfg) {
  prot_dt <- fread(case_cfg$input_path)
  clinical_dt <- fread(case_cfg$clinical_path)
  sample_cols <- names(prot_dt)[-1]
  sample_ids <- sub("^.*:", "", sample_cols)
  keep_sample <- !grepl("^(QC|NCI7-|Pooled Sample$)", sample_ids)
  prot_dt <- prot_dt[, c(1L, which(keep_sample) + 1L), with = FALSE]
  sample_cols <- names(prot_dt)[-1]
  sample_ids <- sub("^.*:", "", sample_cols)
  clinical_dt <- clinical_dt[aliquot_submitter_id %in% sample_ids]
  clinical_dt[, story_group := fifelse(
    tumor_grade %in% case_cfg$grade_low,
    case_cfg$group_levels[1],
    fifelse(tumor_grade %in% case_cfg$grade_high, case_cfg$group_levels[2], NA_character_)
  )]
  clinical_dt <- clinical_dt[!is.na(story_group)]
  keep_ids <- clinical_dt$aliquot_submitter_id
  keep_idx <- sample_ids %in% keep_ids
  prot_dt <- prot_dt[, c(1L, which(keep_idx) + 1L), with = FALSE]
  sample_cols <- names(prot_dt)[-1]
  sample_ids <- sub("^.*:", "", sample_cols)
  clinical_dt <- clinical_dt[match(sample_ids, aliquot_submitter_id)]
  group <- factor(clinical_dt$story_group, levels = case_cfg$group_levels)
  expr <- as.matrix(prot_dt[, -1, with = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- make.unique(prot_dt[[1]])
  expr[is.nan(expr)] <- NA_real_
  n_a <- rowSums(!is.na(expr[, group == levels(group)[1], drop = FALSE]))
  n_b <- rowSums(!is.na(expr[, group == levels(group)[2], drop = FALSE]))
  keep_feature <- n_a >= 2L & n_b >= 2L
  expr <- expr[keep_feature, , drop = FALSE]
  list(
    expr = expr,
    log_expr = expr,
    group = droplevels(group),
    sample_ids = sample_ids,
    clinical = clinical_dt,
    feature_key = "gene",
    assay_name = "proteome"
  )
}

story_read_case_data <- function(case_cfg) {
  switch(
    case_cfg$modality,
    bulk = story_read_bulk_case(case_cfg),
    pseudobulk = story_read_bulk_case(case_cfg),
    proteomics = story_read_proteomics_case(case_cfg)
  )
}

story_case_marker_reference <- function(case_cfg, case_data = NULL) {
  if (!identical(case_cfg$case_id, "pb_kang18")) {
    return(data.table(
      cell_type = character(),
      gene = character(),
      source_label = character(),
      source_url = character()
    ))
  }
  marker_dt <- data.table(
    cell_type = c(
      rep("B cells", 5L),
      rep("CD14+ Monocytes", 5L),
      rep("CD4 T cells", 5L),
      rep("CD8 T cells", 5L),
      rep("Dendritic cells", 5L),
      rep("FCGR3A+ Monocytes", 5L),
      rep("Megakaryocytes", 5L),
      rep("NK cells", 5L)
    ),
    gene = c(
      "MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA",
      "LST1", "S100A8", "S100A9", "FCN1", "CTSS",
      "CD3D", "IL7R", "LTB", "SELL", "GIMAP5",
      "CD8A", "CCL5", "TRAC", "CTSW", "LTB",
      "FCER1A", "CST3", "HLA-DQA1", "GPR183", "HLA-DPA1",
      "FCGR3A", "VMO1", "IFITM3", "SAT1", "IFI30",
      "PPBP", "PF4", "GNG11", "SDPR", "NRGN",
      "GNLY", "NKG7", "GZMB", "FGFBP2", "CTSW"
    ),
    source_label = "Seurat Kang PBMC tutorial canonical markers",
    source_url = "https://satijalab.org/seurat/archive/v3.2/immune_alignment.html"
  )
  if (!is.null(case_data) && !is.null(case_data$cell_types)) {
    marker_dt <- marker_dt[cell_type %in% unique(case_data$cell_types)]
  }
  unique(marker_dt)
}

story_marker_enrichment <- function(standardized_dt, marker_ref_dt) {
  if (nrow(marker_ref_dt) == 0L) {
    return(data.table(
      channel_category = character(),
      cell_type = character(),
      n_channel_genes = integer(),
      n_markers_total = integer(),
      n_markers_in_data = integer(),
      n_overlap = integer(),
      expected_overlap = numeric(),
      fold_enrichment = numeric(),
      p_value = numeric(),
      overlap_genes = character(),
      p.adjust = numeric()
    ))
  }
  universe_genes <- unique(na.omit(standardized_dt$gene))
  categories <- c("DE-only", "DV-only", "DG-only", "overlap")
  out <- rbindlist(lapply(categories, function(cat_name) {
    channel_genes <- unique(standardized_dt[channel_category == cat_name, gene])
    rbindlist(lapply(unique(marker_ref_dt$cell_type), function(cell_name) {
      marker_genes <- unique(marker_ref_dt[cell_type == cell_name, gene])
      marker_genes_in_data <- intersect(marker_genes, universe_genes)
      overlap_genes <- intersect(channel_genes, marker_genes_in_data)
      universe_n <- length(universe_genes)
      channel_n <- length(channel_genes)
      marker_n <- length(marker_genes_in_data)
      overlap_n <- length(overlap_genes)
      expected_overlap <- if (universe_n > 0L) channel_n * marker_n / universe_n else NA_real_
      fold_enrichment <- if (!is.na(expected_overlap) && expected_overlap > 0) overlap_n / expected_overlap else NA_real_
      p_value <- if (universe_n > 0L && channel_n > 0L && marker_n > 0L) {
        stats::phyper(overlap_n - 1L, marker_n, universe_n - marker_n, channel_n, lower.tail = FALSE)
      } else {
        NA_real_
      }
      data.table(
        channel_category = cat_name,
        cell_type = cell_name,
        n_channel_genes = channel_n,
        n_markers_total = length(marker_genes),
        n_markers_in_data = marker_n,
        n_overlap = overlap_n,
        expected_overlap = expected_overlap,
        fold_enrichment = fold_enrichment,
        p_value = p_value,
        overlap_genes = paste(overlap_genes, collapse = ";")
      )
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
  out[, p.adjust := stats::p.adjust(p_value, method = "BH")]
  out[order(channel_category, p.adjust, -n_overlap, cell_type)]
}

story_marker_scores <- function(log_expr, sample_ids, group, marker_ref_dt, clinical = NULL) {
  if (nrow(marker_ref_dt) == 0L) {
    return(data.table(
      sample_id = character(),
      cell_type = character(),
      marker_score = numeric(),
      n_marker_genes = integer(),
      group = character()
    ))
  }
  score_list <- lapply(unique(marker_ref_dt$cell_type), function(cell_name) {
    genes <- intersect(rownames(log_expr), unique(marker_ref_dt[cell_type == cell_name, gene]))
    if (length(genes) == 0L) {
      return(data.table(sample_id = sample_ids, cell_type = cell_name, marker_score = NA_real_, n_marker_genes = 0L))
    }
    mat <- log_expr[genes, , drop = FALSE]
    z_mat <- t(scale(t(mat)))
    z_mat[is.na(z_mat)] <- 0
    data.table(
      sample_id = sample_ids,
      cell_type = cell_name,
      marker_score = as.numeric(colMeans(z_mat)),
      n_marker_genes = length(genes)
    )
  })
  out <- rbindlist(score_list, use.names = TRUE, fill = TRUE)
  out[, group := as.character(group[match(sample_id, sample_ids)])]
  if (!is.null(clinical)) {
    clin_dt <- as.data.table(copy(clinical))
    setnames(clin_dt, old = "aliquot_submitter_id", new = "sample_id", skip_absent = TRUE)
    out <- clin_dt[out, on = "sample_id"]
  }
  out
}

story_representative_marker_map <- function(representative_dt, marker_ref_dt) {
  if (nrow(representative_dt) == 0L || nrow(marker_ref_dt) == 0L) {
    return(data.table(
      gene = character(),
      channel_category = character(),
      rank_within_category = integer(),
      cell_type = character(),
      source_label = character(),
      source_url = character()
    ))
  }
  merge(
    unique(representative_dt[, .(gene, channel_category, rank_within_category)]),
    unique(marker_ref_dt[, .(gene, cell_type, source_label, source_url)]),
    by = "gene",
    all = FALSE,
    allow.cartesian = TRUE
  )
}

story_make_pseudo_counts <- function(log_expr) {
  med_log2 <- median(log_expr, na.rm = TRUE)
  shifted <- log_expr - med_log2 + log2(1000)
  pseudo_counts <- round(2^shifted)
  pseudo_counts[is.na(pseudo_counts) | pseudo_counts < 0] <- 0
  storage.mode(pseudo_counts) <- "integer"
  pseudo_counts
}

story_run_sgcb_case <- function(case_cfg, case_data) {
  if (!is.null(case_cfg$cached_result_path) && isTRUE(case_cfg$use_cached_default)) {
    return(as.data.table(readRDS(case_cfg$cached_result_path)))
  }
  if (case_cfg$modality %in% c("bulk", "pseudobulk")) {
    res <- sgcbDE(counts = case_data$expr, group = case_data$group)
    return(as.data.table(as.data.frame(res)))
  }
  pseudo_counts <- story_make_pseudo_counts(case_data$log_expr)
  res <- sgcbProtein(
    counts = pseudo_counts,
    group = case_data$group,
    prior_de = "auto",
    prior_dv = "auto",
    prior_sd_de = 1.0,
    prior_sd_dv = 0.5
  )
  as.data.table(as.data.frame(res))
}

story_load_or_run_result <- function(case_cfg, case_data, use_cached = TRUE, rerun = FALSE) {
  if (!rerun && use_cached && !is.null(case_cfg$cached_result_path)) {
    return(as.data.table(readRDS(case_cfg$cached_result_path)))
  }
  story_run_sgcb_case(case_cfg, case_data)
}

story_min_with_na <- function(dt, cols) {
  if (length(cols) == 0L) {
    return(rep(NA_real_, nrow(dt)))
  }
  vals <- do.call(pmin, c(dt[, ..cols], list(na.rm = TRUE)))
  vals[is.infinite(vals)] <- NA_real_
  vals
}

story_standardize_result <- function(result_dt, fdr_cutoff = STORY_FDR_CUTOFF) {
  dt <- as.data.table(copy(result_dt))
  if (!"gene" %in% names(dt)) {
    if ("gene_id" %in% names(dt)) {
      dt[, gene := gene_id]
    } else {
      dt[, gene := rownames(result_dt)]
    }
  }
  dt[, gene := as.character(gene)]
  dg_cols <- intersect(c("dg_alpha_padj", "dg_gamma_padj", "dg_shape_padj", "dg_entropy_padj"), names(dt))
  dt[, padj_de := if ("padj" %in% names(dt)) as.numeric(padj) else NA_real_]
  dt[, padj_dv := if ("dv_padj" %in% names(dt)) as.numeric(dv_padj) else if ("padj_dv" %in% names(dt)) as.numeric(padj_dv) else NA_real_]
  dt[, padj_dg := story_min_with_na(dt, dg_cols)]
  dt[, log2FC := if ("log2FoldChange" %in% names(dt)) as.numeric(log2FoldChange) else NA_real_]
  dt[, dv_log2_cv2_ratio := if ("dv_log2_cv2_ratio" %in% names(dt)) as.numeric(dv_log2_cv2_ratio) else NA_real_]
  dt[, dg_alpha_log2_ratio := if ("dg_alpha_log2_ratio" %in% names(dt)) as.numeric(dg_alpha_log2_ratio) else NA_real_]
  dt[, dg_gamma_log2_ratio := if ("dg_gamma_log2_ratio" %in% names(dt)) as.numeric(dg_gamma_log2_ratio) else NA_real_]
  dt[, p_de_post := if ("p_de_post" %in% names(dt)) as.numeric(p_de_post) else NA_real_]
  dt[, p_dv_post := if ("p_dv_post" %in% names(dt)) as.numeric(p_dv_post) else NA_real_]
  dt[, p_dg_post := if ("p_dg_post" %in% names(dt)) as.numeric(p_dg_post) else NA_real_]
  dt[, p_de_only_post := if ("p_de_only_post" %in% names(dt)) as.numeric(p_de_only_post) else NA_real_]
  dt[, p_dv_only_post := if ("p_dv_only_post" %in% names(dt)) as.numeric(p_dv_only_post) else NA_real_]
  dt[, p_dg_only_post := if ("p_dg_only_post" %in% names(dt)) as.numeric(p_dg_only_post) else NA_real_]
  dt[, SGCB_Score := if ("SGCB_Score" %in% names(dt)) as.numeric(SGCB_Score) else NA_real_]
  dt[, sig_de := if ("de_fdr_call" %in% names(dt)) as.logical(de_fdr_call) else (!is.na(padj_de) & padj_de < fdr_cutoff)]
  dt[, sig_dv := if ("dv_fdr_call" %in% names(dt)) as.logical(dv_fdr_call) else (!is.na(padj_dv) & padj_dv < fdr_cutoff)]
  dt[, sig_dg := if ("dg_fdr_call" %in% names(dt)) as.logical(dg_fdr_call) else (!is.na(padj_dg) & padj_dg < fdr_cutoff)]
  dt[is.na(sig_de), sig_de := FALSE]
  dt[is.na(sig_dv), sig_dv := FALSE]
  dt[is.na(sig_dg), sig_dg := FALSE]
  dt[, n_channels := as.integer(sig_de) + as.integer(sig_dv) + as.integer(sig_dg)]
  dt[, channel_signature := fifelse(
    n_channels == 0L,
    "None",
    paste0(
      fifelse(sig_de, "DE", ""),
      fifelse(sig_dv, fifelse(sig_de, "+DV", "DV"), ""),
      fifelse(sig_dg, fifelse(sig_de | sig_dv, "+DG", "DG"), "")
    )
  )]
  dt[, channel_category := fifelse(
    n_channels == 0L,
    "None",
    fifelse(n_channels == 1L, paste0(channel_signature, "-only"), "overlap")
  )]
  dt[, channel_category_detail := fifelse(channel_category == "overlap", channel_signature, channel_category)]
  dt[, de_metric := -log10(pmax(padj_de, 1e-300)) + abs(fcoalesce(log2FC, 0))]
  dt[, dv_metric := -log10(pmax(padj_dv, 1e-300)) + abs(fcoalesce(dv_log2_cv2_ratio, 0))]
  dt[, dg_metric := -log10(pmax(padj_dg, 1e-300)) + abs(fcoalesce(dg_alpha_log2_ratio, 0)) + abs(fcoalesce(dg_gamma_log2_ratio, 0))]
  dt[, overlap_metric := fcoalesce(SGCB_Score, 0)]
  setcolorder(dt, unique(c(
    "gene", "channel_category", "channel_category_detail", "channel_signature",
    "sig_de", "sig_dv", "sig_dg", "padj_de", "padj_dv", "padj_dg",
    "log2FC", "dv_log2_cv2_ratio", "dg_alpha_log2_ratio", "dg_gamma_log2_ratio",
    "p_de_post", "p_dv_post", "p_dg_post", "p_de_only_post", "p_dv_only_post", "p_dg_only_post",
    "SGCB_Score", names(dt)
  )))
  dt
}

story_channel_summary <- function(standardized_dt) {
  standardized_dt[, .(
    n_features = .N,
    median_abs_log2FC = median(abs(log2FC), na.rm = TRUE),
    median_abs_dv_log2_cv2 = median(abs(dv_log2_cv2_ratio), na.rm = TRUE),
    median_abs_dg_gamma_log2 = median(abs(dg_gamma_log2_ratio), na.rm = TRUE),
    median_sgcb_score = median(SGCB_Score, na.rm = TRUE)
  ), by = .(channel_category, channel_category_detail)][order(-n_features, channel_category, channel_category_detail)]
}

story_pick_representative_genes <- function(standardized_dt, top_n = STORY_TOP_FEATURES_PER_CHANNEL) {
  de_dt <- standardized_dt[channel_category == "DE-only"][order(padj_de, -abs(log2FC), gene)]
  dv_dt <- standardized_dt[channel_category == "DV-only"][order(padj_dv, -abs(dv_log2_cv2_ratio), gene)]
  dg_dt <- standardized_dt[channel_category == "DG-only"][order(padj_dg, -abs(dg_gamma_log2_ratio), gene)]
  ov_dt <- standardized_dt[channel_category == "overlap"][order(-overlap_metric, gene)]
  picked <- rbindlist(list(
    de_dt[seq_len(min(top_n, nrow(de_dt))), .(gene, channel_category, rank_metric = de_metric)],
    dv_dt[seq_len(min(top_n, nrow(dv_dt))), .(gene, channel_category, rank_metric = dv_metric)],
    dg_dt[seq_len(min(top_n, nrow(dg_dt))), .(gene, channel_category, rank_metric = dg_metric)],
    ov_dt[seq_len(min(top_n, nrow(ov_dt))), .(gene, channel_category, rank_metric = overlap_metric)]
  ), use.names = TRUE, fill = TRUE)
  picked[, rank_within_category := seq_len(.N), by = channel_category]
  picked
}

story_symbol_mapping <- function(symbols) {
  map <- suppressMessages(clusterProfiler::bitr(unique(symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  as.data.table(unique(map))
}

story_empty_enrichment <- function(category, database, ontology) {
  data.table(
    channel_category = category,
    database = database,
    ontology = ontology,
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    Count = integer(),
    geneID = character(),
    input_n = integer(),
    mapped_n = integer()
  )
}

story_extract_enrichment_dt <- function(enrich_obj, category, database, ontology, input_n, mapped_n) {
  if (is.null(enrich_obj)) {
    out <- story_empty_enrichment(category, database, ontology)
    out[, input_n := input_n]
    out[, mapped_n := mapped_n]
    return(out)
  }
  enrich_df <- as.data.table(as.data.frame(enrich_obj))
  if (nrow(enrich_df) == 0L) {
    out <- story_empty_enrichment(category, database, ontology)
    out[, input_n := input_n]
    out[, mapped_n := mapped_n]
    return(out)
  }
  enrich_df[, channel_category := category]
  enrich_df[, database := database]
  enrich_df[, ontology := ontology]
  enrich_df[, input_n := input_n]
  enrich_df[, mapped_n := mapped_n]
  enrich_df
}

story_run_enrichment <- function(gene_symbols, universe_symbols, category) {
  gene_symbols <- unique(na.omit(gene_symbols))
  universe_symbols <- unique(na.omit(universe_symbols))
  if (length(gene_symbols) < 5L) {
    return(rbindlist(list(
      story_empty_enrichment(category, "GO", "BP"),
      story_empty_enrichment(category, "GO", "MF"),
      story_empty_enrichment(category, "GO", "CC"),
      story_empty_enrichment(category, "KEGG", "KEGG"),
      story_empty_enrichment(category, "Reactome", "Reactome")
    ), use.names = TRUE, fill = TRUE))
  }
  gene_map <- story_symbol_mapping(gene_symbols)
  universe_map <- story_symbol_mapping(universe_symbols)
  gene_entrez <- unique(gene_map$ENTREZID)
  universe_entrez <- unique(universe_map$ENTREZID)
  if (length(gene_entrez) < 5L) {
    return(rbindlist(list(
      story_empty_enrichment(category, "GO", "BP"),
      story_empty_enrichment(category, "GO", "MF"),
      story_empty_enrichment(category, "GO", "CC"),
      story_empty_enrichment(category, "KEGG", "KEGG"),
      story_empty_enrichment(category, "Reactome", "Reactome")
    ), use.names = TRUE, fill = TRUE))
  }
  go_maps <- story_local_go_maps()
  go_bp <- story_local_enricher(
    gene_entrez,
    universe_entrez,
    go_maps$term2gene[ONTOLOGY == "BP", .(GO, ENTREZID)],
    go_maps$term2name[, .(GOID, TERM)]
  )
  go_mf <- story_local_enricher(
    gene_entrez,
    universe_entrez,
    go_maps$term2gene[ONTOLOGY == "MF", .(GO, ENTREZID)],
    go_maps$term2name[, .(GOID, TERM)]
  )
  go_cc <- story_local_enricher(
    gene_entrez,
    universe_entrez,
    go_maps$term2gene[ONTOLOGY == "CC", .(GO, ENTREZID)],
    go_maps$term2name[, .(GOID, TERM)]
  )
  kegg <- tryCatch({
    kegg_maps <- story_local_kegg_maps()
    story_local_enricher(
      gene_entrez,
      universe_entrez,
      kegg_maps$term2gene[, .(pathway_id, ENTREZID)],
      kegg_maps$term2name[, .(pathway_id, pathway_name)]
    )
  }, error = function(e) NULL)
  reactome <- tryCatch(
    ReactomePA::enrichPathway(gene = gene_entrez, universe = universe_entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) NULL
  )
  rbindlist(list(
    story_extract_enrichment_dt(go_bp, category, "GO", "BP", length(gene_symbols), length(gene_entrez)),
    story_extract_enrichment_dt(go_mf, category, "GO", "MF", length(gene_symbols), length(gene_entrez)),
    story_extract_enrichment_dt(go_cc, category, "GO", "CC", length(gene_symbols), length(gene_entrez)),
    story_extract_enrichment_dt(kegg, category, "KEGG", "KEGG", length(gene_symbols), length(gene_entrez)),
    story_extract_enrichment_dt(reactome, category, "Reactome", "Reactome", length(gene_symbols), length(gene_entrez))
  ), use.names = TRUE, fill = TRUE)
}

story_run_case_enrichment <- function(standardized_dt) {
  categories <- c("DE-only", "DV-only", "DG-only", "overlap")
  universe <- standardized_dt$gene
  rbindlist(lapply(categories, function(cat_name) {
    story_run_enrichment(standardized_dt[channel_category == cat_name, gene], universe, cat_name)
  }), use.names = TRUE, fill = TRUE)
}

story_rank_for_signature <- function(dt, category, top_k) {
  ranked <- switch(
    category,
    "DE-only" = dt[channel_category == category][order(padj_de, -abs(log2FC), gene)],
    "DV-only" = dt[channel_category == category][order(padj_dv, -abs(dv_log2_cv2_ratio), gene)],
    "DG-only" = dt[channel_category == category][order(padj_dg, -abs(dg_gamma_log2_ratio), gene)],
    dt[channel_category == category][order(-overlap_metric, gene)]
  )
  ranked[seq_len(min(top_k, nrow(ranked))), gene]
}

story_compute_sample_scores <- function(log_expr, standardized_dt, sample_ids, group, clinical = NULL, top_k = STORY_SIGNATURE_TOP_K) {
  categories <- c("DE-only", "DV-only", "DG-only", "overlap")
  score_list <- lapply(categories, function(cat_name) {
    genes <- intersect(rownames(log_expr), story_rank_for_signature(standardized_dt, cat_name, top_k))
    if (length(genes) == 0L) {
      return(data.table(sample_id = sample_ids, channel_category = cat_name, signature_score = NA_real_, n_signature_genes = 0L))
    }
    mat <- log_expr[genes, , drop = FALSE]
    z_mat <- t(scale(t(mat)))
    z_mat[is.na(z_mat)] <- 0
    scores <- colMeans(z_mat)
    data.table(sample_id = sample_ids, channel_category = cat_name, signature_score = as.numeric(scores), n_signature_genes = length(genes))
  })
  out <- rbindlist(score_list, use.names = TRUE, fill = TRUE)
  out[, group := as.character(group[match(sample_id, sample_ids)])]
  if (!is.null(clinical)) {
    clin_dt <- as.data.table(copy(clinical))
    setnames(clin_dt, old = "aliquot_submitter_id", new = "sample_id", skip_absent = TRUE)
    out <- clin_dt[out, on = "sample_id"]
  }
  out
}

story_prepare_distribution_data <- function(log_expr, sample_ids, group, representative_dt, clinical = NULL) {
  keep_genes <- intersect(rownames(log_expr), unique(representative_dt$gene))
  if (length(keep_genes) == 0L) {
    return(data.table())
  }
  expr_dt <- as.data.table(log_expr[keep_genes, , drop = FALSE], keep.rownames = "gene")
  long_dt <- melt(expr_dt, id.vars = "gene", variable.name = "sample_id", value.name = "expression")
  rep_dt <- unique(representative_dt[, .(gene, channel_category, rank_within_category)])
  long_dt <- rep_dt[long_dt, on = "gene"]
  long_dt[, group := as.character(group[match(sample_id, sample_ids)])]
  if (!is.null(clinical)) {
    clin_dt <- as.data.table(copy(clinical))
    setnames(clin_dt, old = "aliquot_submitter_id", new = "sample_id", skip_absent = TRUE)
    long_dt <- clin_dt[long_dt, on = "sample_id"]
  }
  long_dt[, gene_label := paste0(channel_category, " | ", gene)]
  long_dt
}
