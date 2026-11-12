ROOT <- normalizePath(Sys.getenv("SGCB_REVISION_ROOT"), winslash = "/", mustWork = TRUE)
RECORD_ID <- Sys.getenv("SGCB_STUDY_RECORD_ID","0.5.9")
DATA_ROOT <- normalizePath(Sys.getenv("SGCB_DATA_ROOT",ROOT),winslash="/",mustWork=TRUE)
PREFIX <- paste0("sgcb_",RECORD_ID,"_kirc")
DATA <- file.path(DATA_ROOT, "data")
RESULTS <- file.path(ROOT, "results")
THREADS <- as.integer(Sys.getenv("SGCB_THREADS",as.character(parallel::detectCores())))
Sys.setenv(OMP_NUM_THREADS = THREADS, OPENBLAS_NUM_THREADS = THREADS, MKL_NUM_THREADS = THREADS)
data.table::setDTthreads(THREADS)
.libPaths(c(Sys.getenv("SGCB_LIBRARY",file.path(dirname(ROOT),"library")),.libPaths()))
library(SGCB)
library(data.table)
library(survival)
source(file.path(ROOT,"development","00_current_setup.R"),encoding="UTF-8")

fit <- readRDS(file.path(RESULTS, paste0(PREFIX,"_validation_models.rds")))
cohort <- fread(file.path(RESULTS, paste0(PREFIX,"_patient_scores.csv")))
primary <- readRDS(file.path(RESULTS,paste0(PREFIX,".rds")))
verify_study_record(file.path(RESULTS,paste0(PREFIX,".rds")),attr(primary,"native_source_id"),RECORD_ID)
stopifnot(identical(fit$run_id,PREFIX),identical(fit$source_id,attr(primary,"native_source_id")))
discovery <- as.data.table(primary)
ds <- readRDS(file.path(DATA, "tcga_kirc_cache", "tcga_kirc.rds"))
expr <- ds$expression[, substr(colnames(ds$expression), 14, 15) == "01", drop = FALSE]
patient <- substr(colnames(expr), 1, 12)
logx <- log2(sweep(expr, 2, colSums(expr), "/") * 1e6 + 0.5)
parameters <- fit$score_parameters
genes <- rownames(parameters$weights)[parameters$weights[, "DV_only"] > 0]
z <- (logx[genes, , drop = FALSE] - parameters$gene_mean[genes]) / parameters$gene_sd[genes]
direction <- sign(discovery[match(genes, gene_id), log2FoldChange])
signed <- colMeans(z * direction)
absolute <- colMeans(abs(z))
stopifnot(identical(cohort$patient[match(patient, cohort$patient)], patient))
training <- cohort$training[match(patient, cohort$patient)]
alternatives <- data.table(patient = patient,
    signed_mean = signed / sd(signed[training]),
    absolute_deviation = absolute / sd(absolute[training]))
cohort <- merge(cohort, alternatives, by = "patient")
cohort[, `:=`(stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
              grade = factor(grade, levels = c("G1/G2", "G3/G4")))]
validation <- cohort[validation == TRUE & is.finite(time) & time > 0 & !is.na(event) &
                       complete.cases(cohort[, .(age, stage, grade)])]
model_signed <- coxph(Surv(time, event) ~ signed_mean + age + strata(stage) + grade,
                      data = validation, x = TRUE)
model_absolute <- coxph(Surv(time, event) ~ absolute_deviation + age + strata(stage) + grade,
                        data = validation, x = TRUE)
table <- rbindlist(list(
  as.data.table(broom::tidy(model_signed, conf.int = TRUE, exponentiate = TRUE))[term == "signed_mean"],
  as.data.table(broom::tidy(model_absolute, conf.int = TRUE, exponentiate = TRUE))[term == "absolute_deviation"]
))
table[, `:=`(qvalue = p.adjust(p.value, "BH"), n = nrow(validation), events = sum(validation$event))]
fwrite(table, file.path(RESULTS, paste0(PREFIX,"_alternative_scores_survival.csv")))
fwrite(alternatives, file.path(RESULTS, paste0(PREFIX,"_alternative_scores.csv")))
fwrite(rbindlist(list(
  as.data.table(cox.zph(model_signed)$table, keep.rownames = "term")[, score := "signed_mean"],
  as.data.table(cox.zph(model_absolute)$table, keep.rownames = "term")[, score := "absolute_deviation"]
)), file.path(RESULTS, paste0(PREFIX,"_alternative_scores_PH.csv")))

external <- readRDS(file.path(RESULTS, paste0(PREFIX,"_external_validation.rds")))
stopifnot(identical(external$run_id,PREFIX),identical(external$source_id,attr(primary,"native_source_id")))
external_cohort <- fread(file.path(RESULTS, paste0(PREFIX,"_external_scores.csv")))
external_genes <- external$signatures[channel == "DV_only", gene_id]
external_x <- external$expression[external_genes, , drop = FALSE]
external_z <- (external_x - rowMeans(external_x)) / matrixStats::rowSds(external_x)
external_direction <- sign(discovery[match(external_genes, gene_id), log2FoldChange])
external_signed <- colMeans(external_z * external_direction)
external_absolute <- colMeans(abs(external_z))
external_scores <- data.table(patient = colnames(external_x),
    signed_mean = external_signed / sd(external_signed),
    absolute_deviation = external_absolute / sd(external_absolute))
external_cohort <- merge(external_cohort, external_scores, by = "patient")
external_complete <- external_cohort[complete.cases(external_cohort[, .(age, stage_high, grade_high)])]
external_signed_model <- coxph(Surv(time, event) ~ signed_mean + age + stage_high + grade_high,
                               data = external_complete, x = TRUE)
external_absolute_model <- coxph(Surv(time, event) ~ absolute_deviation + age + stage_high + grade_high,
                                 data = external_complete, x = TRUE)
external_table <- rbindlist(list(
  as.data.table(broom::tidy(external_signed_model, conf.int = TRUE, exponentiate = TRUE))[term == "signed_mean"],
  as.data.table(broom::tidy(external_absolute_model, conf.int = TRUE, exponentiate = TRUE))[term == "absolute_deviation"]
))
external_table[, `:=`(qvalue = p.adjust(p.value, "BH"), n = nrow(external_complete),
                      events = sum(external_complete$event))]
fwrite(external_table, file.path(RESULTS, paste0(PREFIX,"_external_alternative_scores_survival.csv")))
fwrite(external_scores, file.path(RESULTS, paste0(PREFIX,"_external_alternative_scores.csv")))
fwrite(rbindlist(list(
  as.data.table(cox.zph(external_signed_model)$table, keep.rownames = "term")[, score := "signed_mean"],
  as.data.table(cox.zph(external_absolute_model)$table, keep.rownames = "term")[, score := "absolute_deviation"]
)), file.path(RESULTS, paste0(PREFIX,"_external_alternative_scores_PH.csv")))

# Figure examples are chosen only by abundance within the three discovery sets.
example_genes <- discovery[channel != "neither"][order(channel, -baseMean, gene_id), .SD[1], by = channel]
input <- readRDS(file.path(DATA, "kirc_independent_discovery.rds"))
x <- input$counts[match(discovery$gene_id, rownames(input$counts)), , drop = FALSE]
libraries <- colSums(x)
tmm <- SGCB:::tmm_factors_cpp(x, as.integer(which.min(abs(libraries - median(libraries))) - 1L),
                             0.3, 0.05, 1e-8, 0L)
effective <- libraries * tmm
sf <- effective / exp(mean(log(effective)))
normalized <- sweep(x, 2, sf, "/") + 0.5
examples <- as.data.table(log2(normalized[match(example_genes$gene_id, discovery$gene_id), , drop = FALSE]),
                         keep.rownames = "gene_id")
examples[, gene_id := example_genes$gene_id]
examples <- melt(examples, id.vars = "gene_id", variable.name = "sample", value.name = "log2_abundance")
examples[, group := input$group[match(sample, colnames(input$counts))]]
examples <- merge(examples, example_genes[, .(gene_id, channel, baseMean, padj, dv_padj)], by = "gene_id")
fwrite(examples, file.path(RESULTS, paste0(PREFIX,"_figure_examples.csv")))
print(table[, .(term, estimate, conf.low, conf.high, p.value, qvalue)])
print(external_table[, .(term, estimate, conf.low, conf.high, p.value, qvalue)])
print(example_genes[, .(gene_id, channel, baseMean, padj, dv_padj)])
