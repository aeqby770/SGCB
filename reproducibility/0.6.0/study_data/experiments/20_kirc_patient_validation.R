ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"), winslash = "/", mustWork = TRUE)
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.5.9")
PREFIX <- paste0("sgcb_",RECORD_ID,"_kirc")
DATA <- file.path(DATA_ROOT, "data")
RESULTS <- file.path(ROOT, "results")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
IMMUNE_GENES <- c("PTPRC", "CD3D", "CD3E", "CD3G", "CD8A", "LST1", "TYROBP")
CHANNELS <- c("DE_only", "DV_only", "DE_and_DV")
Sys.setenv(OMP_NUM_THREADS = THREADS, OPENBLAS_NUM_THREADS = THREADS, MKL_NUM_THREADS = THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
library(survival)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")

primary <- readRDS(file.path(RESULTS,paste0(PREFIX,".rds")))
stopifnot(identical(attr(primary,"native_source_id"),SOURCE_ID))
res <- as.data.table(primary)
membership <- fread(file.path(DATA_ROOT, "results", "kirc_discovery_sample_membership.csv"))
ds <- readRDS(file.path(DATA, "tcga_kirc_cache", "tcga_kirc.rds"))
clinical <- as.data.table(readRDS(file.path(DATA, "tcga_kirc_cache", "tcga_kirc_clinical.rds")))
expr <- ds$expression[, substr(colnames(ds$expression), 14, 15) == "01", drop = FALSE]
pid <- substr(colnames(expr), 1, 12)
stopifnot(!anyDuplicated(pid))
train <- pid %in% membership[group == "tumor", patient]
validation <- !pid %in% membership$patient
logx <- log2(sweep(expr, 2, colSums(expr), "/") * 1e6 + 0.5)
train_mean <- rowMeans(logx[, train, drop = FALSE])
train_sd <- matrixStats::rowSds(logx[, train, drop = FALSE])
genes <- intersect(res[channel %in% CHANNELS, gene_id], rownames(logx)[train_sd > 0])
strata <- res[match(genes, gene_id), channel]
weights <- 1.0 * outer(strata, CHANNELS, "==")
weights <- sweep(weights, 2, colSums(weights), "/")
rownames(weights) <- genes
colnames(weights) <- CHANNELS
z <- (logx[genes, , drop = FALSE] - train_mean[genes]) / train_sd[genes]
score <- t(z) %*% weights
score_sd <- matrixStats::colSds(score[train, , drop = FALSE])
score_cut <- matrixStats::colMedians(score[train, , drop = FALSE])
standard_score <- sweep(score, 2, score_sd, "/")
immune_z <- (logx[IMMUNE_GENES, , drop = FALSE] - train_mean[IMMUNE_GENES]) / train_sd[IMMUNE_GENES]
immune <- colMeans(immune_z)
immune <- immune / sd(immune[train])
published <- fread(file.path(DATA, "kirc_published_purity.csv"))
published[, patient := substr(`Sample ID`, 1, 12)]
purity <- published[, .(CPE = median(CPE, na.rm = TRUE)), by = patient]
scores <- data.table(patient = pid, sample = colnames(expr), training = train,
                     validation = validation, standard_score, immune_proxy = immune,
                     DV_high = factor(fifelse(score[, "DV_only"] > score_cut[2], "High", "Low"),
                                      levels = c("Low", "High")))
clinical[, `:=`(patient = submitter_id,
    time = fifelse(vital_status == "Dead", as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
    event = as.integer(vital_status == "Dead"), age = as.numeric(age_at_diagnosis) / 365.25,
    stage = factor(ajcc_pathologic_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    grade = factor(fifelse(tumor_grade %in% c("G1", "G2"), "G1/G2",
                   fifelse(tumor_grade %in% c("G3", "G4"), "G3/G4", NA_character_))))]
cohort <- merge(scores, clinical[, .(patient, time, event, age, stage, grade)], by = "patient", all.x = TRUE)
cohort <- merge(cohort, purity, by = "patient", all.x = TRUE)
validation_data <- cohort[validation == TRUE & is.finite(time) & time > 0 & !is.na(event)]
complete <- validation_data[complete.cases(validation_data[, .(age, stage, grade)])]
composition <- complete[is.finite(CPE) & is.finite(immune_proxy)]

de_u <- coxph(Surv(time, event) ~ DE_only, data = validation_data, x = TRUE)
dv_u <- coxph(Surv(time, event) ~ DV_only, data = validation_data, x = TRUE)
both_u <- coxph(Surv(time, event) ~ DE_and_DV, data = validation_data, x = TRUE)
de_a <- coxph(Surv(time, event) ~ DE_only + age + strata(stage) + grade, data = complete, x = TRUE)
dv_a <- coxph(Surv(time, event) ~ DV_only + age + strata(stage) + grade, data = complete, x = TRUE)
both_a <- coxph(Surv(time, event) ~ DE_and_DV + age + strata(stage) + grade, data = complete, x = TRUE)
de_c <- coxph(Surv(time, event) ~ DE_only + age + strata(stage) + grade + CPE + immune_proxy,
              data = composition, x = TRUE)
dv_c <- coxph(Surv(time, event) ~ DV_only + age + strata(stage) + grade + CPE + immune_proxy,
              data = composition, x = TRUE)
both_c <- coxph(Surv(time, event) ~ DE_and_DV + age + strata(stage) + grade + CPE + immune_proxy,
                data = composition, x = TRUE)
tab <- rbindlist(list(
  as.data.table(broom::tidy(de_u, exponentiate = TRUE, conf.int = TRUE))[term == "DE_only"][, model := "unadjusted"],
  as.data.table(broom::tidy(dv_u, exponentiate = TRUE, conf.int = TRUE))[term == "DV_only"][, model := "unadjusted"],
  as.data.table(broom::tidy(both_u, exponentiate = TRUE, conf.int = TRUE))[term == "DE_and_DV"][, model := "unadjusted"],
  as.data.table(broom::tidy(de_a, exponentiate = TRUE, conf.int = TRUE))[term == "DE_only"][, model := "clinical"],
  as.data.table(broom::tidy(dv_a, exponentiate = TRUE, conf.int = TRUE))[term == "DV_only"][, model := "clinical"],
  as.data.table(broom::tidy(both_a, exponentiate = TRUE, conf.int = TRUE))[term == "DE_and_DV"][, model := "clinical"],
  as.data.table(broom::tidy(de_c, exponentiate = TRUE, conf.int = TRUE))[term == "DE_only"][, model := "composition"],
  as.data.table(broom::tidy(dv_c, exponentiate = TRUE, conf.int = TRUE))[term == "DV_only"][, model := "composition"],
  as.data.table(broom::tidy(both_c, exponentiate = TRUE, conf.int = TRUE))[term == "DE_and_DV"][, model := "composition"]
))
tab[, qvalue := p.adjust(p.value, "BH"), by = model]
size <- data.table(model = c("unadjusted", "clinical", "composition"),
    n = c(nrow(validation_data), nrow(complete), nrow(composition)),
    events = c(sum(validation_data$event), sum(complete$event), sum(composition$event)))
tab <- merge(tab, size, by = "model")
fwrite(tab, file.path(RESULTS, paste0(PREFIX,"_heldout_survival.csv")))
ph <- rbindlist(list(
  as.data.table(cox.zph(de_a)$table, keep.rownames = "term")[, signature := "DE_only"],
  as.data.table(cox.zph(dv_a)$table, keep.rownames = "term")[, signature := "DV_only"],
  as.data.table(cox.zph(both_a)$table, keep.rownames = "term")[, signature := "DE_and_DV"]
))
fwrite(ph, file.path(RESULTS, paste0(PREFIX,"_heldout_proportional_hazards.csv")))
km <- survfit(Surv(time, event) ~ DV_high, data = validation_data)
fwrite(as.data.table(broom::tidy(km)), file.path(RESULTS, paste0(PREFIX,"_heldout_KM.csv")))
km_cox <- coxph(Surv(time, event) ~ DV_high, data = validation_data)
km_logrank <- survdiff(Surv(time, event) ~ DV_high, data = validation_data)
km_tab <- as.data.table(broom::tidy(km_cox, exponentiate = TRUE, conf.int = TRUE))
km_tab[, `:=`(cutoff = score_cut[2] / score_sd[2], logrank_p = pchisq(km_logrank$chisq, 1, lower.tail = FALSE),
              n_low = sum(validation_data$DV_high == "Low"), n_high = sum(validation_data$DV_high == "High"))]
fwrite(km_tab, file.path(RESULTS, paste0(PREFIX,"_heldout_KM_test.csv")))
fwrite(cohort, file.path(RESULTS, paste0(PREFIX,"_patient_scores.csv")))
fwrite(data.table(gene_id = genes, channel = strata, training_mean = train_mean[genes],
                  training_sd = train_sd[genes]), file.path(RESULTS, paste0(PREFIX,"_locked_signature_genes.csv")))
fwrite(data.table(channel = CHANNELS, genes = colSums(weights > 0),
                  training_score_sd = score_sd, training_median_cutoff = score_cut),
       file.path(RESULTS, paste0(PREFIX,"_signature_parameters.csv")))
fwrite(data.table(training_tumors = sum(train), discovery_patients_excluded = uniqueN(membership$patient),
                  primary_tumors = ncol(expr), heldout_primary_tumors = sum(validation),
                  heldout_survival = nrow(validation_data), heldout_events = sum(validation_data$event),
                  rho_DV_CPE = cor(composition$DV_only, composition$CPE, method = "spearman"),
                  rho_DV_immune = cor(composition$DV_only, composition$immune_proxy, method = "spearman")),
       file.path(RESULTS, paste0(PREFIX,"_validation_summary.csv")))
saveRDS(list(run_id=PREFIX,source_id=SOURCE_ID,
             models = list(de_a, dv_a, both_a, de_c, dv_c, both_c), km = km,
             score_parameters = list(gene_mean = train_mean, gene_sd = train_sd,
                                     weights = weights, score_sd = score_sd, cutoff = score_cut)),
        file.path(RESULTS, paste0(PREFIX,"_validation_models.rds")))
print(tab[, .(model, term, estimate, conf.low, conf.high, p.value, qvalue, n, events)])
print(km_tab)
