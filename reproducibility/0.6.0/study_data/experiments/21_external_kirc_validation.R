ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"), winslash = "/", mustWork = TRUE)
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.5.9")
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
PREFIX <- paste0("sgcb_",RECORD_ID,"_kirc")
DATA <- file.path(DATA_ROOT, "data")
RESULTS <- file.path(ROOT, "results")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
CHANNELS <- c("DE_only", "DV_only", "DE_and_DV")
Sys.setenv(OMP_NUM_THREADS = THREADS, OPENBLAS_NUM_THREADS = THREADS, MKL_NUM_THREADS = THREADS)
data.table::setDTthreads(THREADS)
library(data.table)
library(survival)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")

raw <- fread(file.path(DATA, "ccRCC_exp_log_quantile_normalized.txt"))
samples <- grep("^ccRCC-", names(raw), value = TRUE)
raw[, refseq := sub("\\..*$", "", SystematicName)]
annotation <- as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(raw$refseq), keytype = "REFSEQ", columns = "SYMBOL"))
annotation <- unique(annotation[!is.na(SYMBOL) & SYMBOL != ""])
annotation <- annotation[, .(SYMBOL = SYMBOL[1], n_symbols = uniqueN(SYMBOL)), by = REFSEQ][n_symbols == 1]
mapped <- merge(raw[, c("REF", "refseq", samples), with = FALSE], annotation,
                 by.x = "refseq", by.y = "REFSEQ", all = FALSE)
expr_sum <- rowsum(as.matrix(mapped[, ..samples]), group = mapped$SYMBOL, reorder = TRUE)
expr <- expr_sum / as.numeric(table(mapped$SYMBOL)[rownames(expr_sum)])
primary <- readRDS(file.path(RESULTS,paste0(PREFIX,".rds")))
verify_study_record(file.path(RESULTS,paste0(PREFIX,".rds")),attr(primary,"native_source_id"),RECORD_ID)
res <- as.data.table(primary)
sd_external <- matrixStats::rowSds(expr)
genes <- intersect(res[channel %in% CHANNELS, gene_id], rownames(expr)[sd_external > 0])
strata <- res[match(genes, gene_id), channel]
weights <- 1.0 * outer(strata, CHANNELS, "==")
weights <- sweep(weights, 2, colSums(weights), "/")
colnames(weights) <- CHANNELS
z <- (expr[genes, , drop = FALSE] - rowMeans(expr[genes, , drop = FALSE])) / sd_external[genes]
# Platform-specific centering uses expression only; gene memberships and
# equal positive weights remain fixed from the TCGA discovery contrast.
scores <- t(z) %*% weights
scores <- sweep(scores, 2, matrixStats::colSds(scores), "/")
clinical <- fread(file.path(DATA, "sato_external_clinical.csv"))
clinical[, `:=`(patient = `sample ID`, time = as.numeric(`observation period (month)`),
                event = as.integer(outcome == "dead"), age = as.numeric(Age),
                T_stage = as.integer(sub("^.*T([1-4]).*$", "\\1", `Stage at diagnosis`)),
                grade = as.numeric(`Fuhrman grade`))]
clinical[, `:=`(stage_high = as.integer(T_stage >= 3), grade_high = as.integer(grade >= 3))]
cohort <- merge(data.table(patient = samples, scores),
    clinical[, .(patient, time, event, age, T_stage, stage_high, grade, grade_high)], by = "patient")
cohort <- cohort[is.finite(time) & time > 0 & !is.na(event)]
complete <- cohort[complete.cases(cohort[, .(age, stage_high, grade_high)])]
de_u <- coxph(Surv(time, event) ~ DE_only, data = cohort, x = TRUE)
dv_u <- coxph(Surv(time, event) ~ DV_only, data = cohort, x = TRUE)
both_u <- coxph(Surv(time, event) ~ DE_and_DV, data = cohort, x = TRUE)
de_a <- coxph(Surv(time, event) ~ DE_only + age + stage_high + grade_high, data = complete, x = TRUE)
dv_a <- coxph(Surv(time, event) ~ DV_only + age + stage_high + grade_high, data = complete, x = TRUE)
both_a <- coxph(Surv(time, event) ~ DE_and_DV + age + stage_high + grade_high, data = complete, x = TRUE)
tab <- rbindlist(list(
  as.data.table(broom::tidy(de_u, exponentiate = TRUE, conf.int = TRUE))[term == "DE_only"][, model := "unadjusted"],
  as.data.table(broom::tidy(dv_u, exponentiate = TRUE, conf.int = TRUE))[term == "DV_only"][, model := "unadjusted"],
  as.data.table(broom::tidy(both_u, exponentiate = TRUE, conf.int = TRUE))[term == "DE_and_DV"][, model := "unadjusted"],
  as.data.table(broom::tidy(de_a, exponentiate = TRUE, conf.int = TRUE))[term == "DE_only"][, model := "clinical"],
  as.data.table(broom::tidy(dv_a, exponentiate = TRUE, conf.int = TRUE))[term == "DV_only"][, model := "clinical"],
  as.data.table(broom::tidy(both_a, exponentiate = TRUE, conf.int = TRUE))[term == "DE_and_DV"][, model := "clinical"]
))
tab[, qvalue := p.adjust(p.value, "BH"), by = model]
tab[, `:=`(n = fifelse(model == "clinical", nrow(complete), nrow(cohort)),
            events = fifelse(model == "clinical", sum(complete$event), sum(cohort$event)))]
fwrite(tab, file.path(RESULTS, paste0(PREFIX,"_external_survival.csv")))
fwrite(cohort, file.path(RESULTS, paste0(PREFIX,"_external_scores.csv")))
fwrite(data.table(channel = CHANNELS, mapped_genes = colSums(weights > 0),
    discovery_genes = as.integer(table(factor(res$channel, levels = CHANNELS))),
    n_expression_samples = ncol(expr), n_survival = nrow(cohort), events = sum(cohort$event),
    annotation_version = as.character(packageVersion("org.Hs.eg.db"))),
    file.path(RESULTS, paste0(PREFIX,"_external_coverage.csv")))
ph <- rbindlist(list(
  as.data.table(cox.zph(de_a)$table, keep.rownames = "term")[, signature := "DE_only"],
  as.data.table(cox.zph(dv_a)$table, keep.rownames = "term")[, signature := "DV_only"],
  as.data.table(cox.zph(both_a)$table, keep.rownames = "term")[, signature := "DE_and_DV"]
))
fwrite(ph, file.path(RESULTS, paste0(PREFIX,"_external_proportional_hazards.csv")))
saveRDS(list(run_id=PREFIX,source_id=attr(primary,"native_source_id"),
    analysis_source_id=SOURCE_ID,analysis_package_version=as.character(packageVersion("SGCB")),
             expression = expr, mapping = annotation, signatures = data.table(gene_id = genes, channel = strata),
             models = list(de_u, dv_u, both_u, de_a, dv_a, both_a)),
        file.path(RESULTS, paste0(PREFIX,"_external_validation.rds")))
print(tab[, .(model, term, estimate, conf.low, conf.high, p.value, qvalue, n, events)])
