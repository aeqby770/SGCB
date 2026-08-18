# =============================================================================
# =============================================================================

loaded <- switch(DATASET,

  # ----- CPTAC Study 6 -----
  "cptac" = {
    pep_dt <- fread(file.path(DATA_DIR, "cptac_lab3_peptides.txt"), integer64 = "numeric")
    pep_dt <- pep_dt[is.na(Reverse) | Reverse != "+"]
    pep_dt <- pep_dt[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]
    int_cols <- grep("^Intensity ", names(pep_dt), value = TRUE)
    int_mat  <- as.matrix(pep_dt[, ..int_cols])
    int_mat[int_mat == 0] <- NA
    prot_id <- pep_dt[["Leading razor protein"]]
    log2_pep <- log2(int_mat)
    pep_tbl  <- data.table(protein = prot_id, log2_pep)
    pep_tbl  <- pep_tbl[!is.na(protein) & protein != ""]
    pep_counts <- pep_tbl[, .(n_pep = .N), by = protein]
    prot_long <- melt(pep_tbl, id.vars = "protein", variable.name = "sample", value.name = "log2int")
    prot_wide <- dcast(prot_long, protein ~ sample, value.var = "log2int", fun.aggregate = median, na.rm = TRUE)
    prot_mat  <- as.matrix(prot_wide[, -1])
    rownames(prot_mat) <- prot_wide$protein
    prot_mat[is.nan(prot_mat)] <- NA
    grp_vec   <- fifelse(grepl("6A", colnames(prot_mat)), "6A", "6B")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "6A", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "6B", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    is_spike <- grepl("HUMAN|UPS", rownames(prot_mat))
    pep_count_matched <- pep_counts[match(rownames(prot_mat), pep_counts$protein), n_pep]
    pep_count_matched[is.na(pep_count_matched)] <- 1L
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("6A", "6B")),
         pep_count = pep_count_matched, is_spike = is_spike)
  },

  # ----- PXD000279 E.coli LFQ -----
  "ecoli_lfq" = {
    pg <- fread(file = file.path(DATA_DIR, "ecoli_lfq_deqms/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(Contaminant) | Contaminant != "+"]
    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    grp_vec <- fifelse(grepl("^LFQ intensity H", lfq_cols), "H", "L")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "H", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "L", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    fasta <- pg[["Fasta headers"]][keep]
    is_spike <- grepl("ECOLI", fasta, ignore.case = TRUE)
    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("L", "H")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- PXD002370 Ratio2 -----
  "sgsds_ratio2" = {
    pg <- fread(file = file.path(DATA_DIR, "sgsds_pxd002370/txt/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]
    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    grp_vec <- fifelse(grepl("^LFQ intensity D", lfq_cols), "D", "E")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "D", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "E", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    fasta <- pg[["Fasta headers"]][keep]
    is_spike <- grepl("HUMAN|UPS", fasta, ignore.case = TRUE)
    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("D", "E")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- PXD002370 Ratio2.5 -----
  "sgsds_ratio2.5" = {
    pg <- fread(file = file.path(DATA_DIR, "sgsds_pxd002370/ratio2.5/txt/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]
    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    grp_vec <- fifelse(grepl("^LFQ intensity C", lfq_cols), "C", "D")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "C", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "D", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    fasta <- pg[["Fasta headers"]][keep]
    is_spike <- grepl("HUMAN|UPS", fasta, ignore.case = TRUE)
    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("C", "D")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- CPTAC UPS1 protein-level (6A vs 6B, UPS spike-in) -----
  "cptac_ups1" = {
    pg <- fread(file = file.path(DATA_DIR, "cptac_ups1/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]
    lfq_cols <- grep("^LFQ intensity ", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority protein IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    grp_vec <- fifelse(grepl("6A", lfq_cols), "6A", "6B")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "6A", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "6B", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    is_spike <- grepl("UPS", rownames(prot_mat))
    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("6A", "6B")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- O'Brien 2018 3-species (condition1 vs condition2, YEAST+ECOLI spike-in) -----
  "obrien_3species" = {
    pg <- fread(file = file.path(DATA_DIR, "obrien2018_3species/proteinGroups.txt"), integer64 = "numeric")
    pg <- pg[is.na(Reverse) | Reverse != "+"]
    pg <- pg[is.na(Potential.contaminant) | Potential.contaminant != "+"]
    lfq_cols <- grep("^LFQ\\.intensity\\.", names(pg), value = TRUE)
    prot_mat <- as.matrix(pg[, ..lfq_cols])
    rownames(prot_mat) <- pg[["Majority.protein.IDs"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    # 2 conditions: first 6 samples (CG1407-CG51963) vs last 6 (GsbPI-S2R) based on README
    grp_vec <- fifelse(grepl("CG1407|CG4676|CG51963|CG5620A|CG5620B|CG5880", lfq_cols), "cond1", "cond2")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "cond1", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "cond2", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    fasta <- pg[["Fasta.headers"]][keep]
    is_spike <- grepl("YEAST|ECOLI", fasta, ignore.case = TRUE)
    pep_cnt <- pg[["Peptides"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("cond1", "cond2")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  },

  # ----- Simulated proteomics -----
  "sim_proteomics" = {
    sim_d <- readRDS(file.path(DATA_DIR, "sim_proteomics.rds"))
    list(prot_log2 = sim_d$prot_log2, group = sim_d$group,
         pep_count = sim_d$pep_count, is_spike = sim_d$is_spike,
         sim_truth = sim_d$truth)
  },

  # ----- PXD004163 TMT miR -----
  "tmt_mir" = {
    pg <- fread(file = file.path(DATA_DIR, "tmt_mir_deqms/Yan_miR_Protein_table.txt"), integer64 = "numeric")
    tmt_cols <- grep("tmt10plex_1[23]", names(pg), value = TRUE)
    tmt_cols <- tmt_cols[!grepl("quanted", tmt_cols)]
    prot_mat <- as.matrix(pg[, ..tmt_cols])
    rownames(prot_mat) <- pg[["Protein accession"]]
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    grp_vec <- fifelse(grepl("126|127N|127C|128N", tmt_cols), "ctrl", "miR372")
    n_va <- rowSums(!is.na(prot_mat[, grp_vec == "ctrl", drop = FALSE]))
    n_vb <- rowSums(!is.na(prot_mat[, grp_vec == "miR372", drop = FALSE]))
    keep <- n_va >= 2 & n_vb >= 2
    prot_mat <- prot_mat[keep, ]
    pep_cnt <- pg[["miR FASP_# PSMs"]][keep]
    pep_cnt[is.na(pep_cnt)] <- 1L
    is_spike <- rep(FALSE, sum(keep))
    list(prot_log2 = prot_mat, group = factor(grp_vec, levels = c("ctrl", "miR372")),
         pep_count = as.integer(pep_cnt), is_spike = is_spike)
  }
)

cat(sprintf("  Loaded %s: %d proteins x %d samples, spike-in: %d\n",
            DATASET, nrow(loaded$prot_log2), ncol(loaded$prot_log2), sum(loaded$is_spike)))
