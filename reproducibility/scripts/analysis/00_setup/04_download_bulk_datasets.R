# =============================================================================
# =============================================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)

# =============================================================================
# =============================================================================

DATA_DIR <- paste0(PROJECT_ROOT, "/benchmark/data")
SC_DIR <- file.path(DATA_DIR, "single_cell")
BULK_DIR <- file.path(DATA_DIR, "bulk")

dir.create(SC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(BULK_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# =============================================================================

cat("\n========== ==========\n")

BiocManager::install(c(
    "scRNAseq",
    "TENxPBMCData
", 
    "muscData",
    "celldex",
    "scater",
    "SingleCellExperiment",
    "ExperimentHub",
    "AnnotationHub",
    "recount3",
    "curatedTCGAData",
    "TCGAutils",
    "GEOquery",
    "SummarizedExperiment",
    "MultiAssayExperiment",
    "BiocFileCache"
), update = FALSE, ask = FALSE)

library(scRNAseq)
library(TENxPBMCData)
library(muscData)
library(recount3)
library(curatedTCGAData)
library(TCGAutils)
library(GEOquery)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentHub)

# =============================================================================
# =============================================================================

cat("\n========== : scRNAseq ==========\n")

cat("\n--- Baron Pancreas (Human) ---\n")
baron_human <- BaronPancreasData("human")
saveRDS(baron_human, file.path(SC_DIR, "baron_pancreas_human.rds"))
cat("Saved: ", ncol(baron_human), " cells, ", nrow(baron_human), " genes\n")

cat("\n--- Baron Pancreas (Mouse) ---\n")
baron_mouse <- BaronPancreasData("mouse")
saveRDS(baron_mouse, file.path(SC_DIR, "baron_pancreas_mouse.rds"))
cat("Saved: ", ncol(baron_mouse), " cells, ", nrow(baron_mouse), " genes\n")

cat("\n--- Segerstolpe Pancreas ---\n")
segerstolpe <- SegerstolpePancreasData()
saveRDS(segerstolpe, file.path(SC_DIR, "segerstolpe_pancreas.rds"))
cat("Saved: ", ncol(segerstolpe), " cells, ", nrow(segerstolpe), " genes\n")

cat("\n--- Muraro Pancreas ---\n")
muraro <- MuraroPancreasData()
saveRDS(muraro, file.path(SC_DIR, "muraro_pancreas.rds"))
cat("Saved: ", ncol(muraro), " cells, ", nrow(muraro), " genes\n")

cat("\n--- Lawlor Pancreas ---\n")
lawlor <- LawlorPancreasData()
saveRDS(lawlor, file.path(SC_DIR, "lawlor_pancreas.rds"))
cat("Saved: ", ncol(lawlor), " cells, ", nrow(lawlor), " genes\n")

cat("\n--- Zeisel Brain ---\n")
zeisel <- ZeiselBrainData()
saveRDS(zeisel, file.path(SC_DIR, "zeisel_brain.rds"))
cat("Saved: ", ncol(zeisel), " cells, ", nrow(zeisel), " genes\n")

cat("\n--- Tasic Brain ---\n")
tasic <- TasicBrainData()
saveRDS(tasic, file.path(SC_DIR, "tasic_brain.rds"))
cat("Saved: ", ncol(tasic), " cells, ", nrow(tasic), " genes\n")

cat("\n--- Campbell Brain ---\n")
campbell <- CampbellBrainData()
saveRDS(campbell, file.path(SC_DIR, "campbell_brain.rds"))
cat("Saved: ", ncol(campbell), " cells, ", nrow(campbell), " genes\n")

cat("\n--- Romanov Brain ---\n")
romanov <- RomanovBrainData()
saveRDS(romanov, file.path(SC_DIR, "romanov_brain.rds"))
cat("Saved: ", ncol(romanov), " cells, ", nrow(romanov), " genes\n")

cat("\n--- MacParland Liver ---\n")
macparland <- MacParlandLiverData()
saveRDS(macparland, file.path(SC_DIR, "macparland_liver.rds"))
cat("Saved: ", ncol(macparland), " cells, ", nrow(macparland), " genes\n")

cat("\n--- Grun Pancreas ---\n")
grun <- GrunPancreasData()
saveRDS(grun, file.path(SC_DIR, "grun_pancreas.rds"))
cat("Saved: ", ncol(grun), " cells, ", nrow(grun), " genes\n")

cat("\n--- Macosko Retina ---\n")
macosko <- MacoskoRetinaData()
saveRDS(macosko, file.path(SC_DIR, "macosko_retina.rds"))
cat("Saved: ", ncol(macosko), " cells, ", nrow(macosko), " genes\n")

cat("\n--- Shekhar Retina ---\n")
shekhar <- ShekharRetinaData()
saveRDS(shekhar, file.path(SC_DIR, "shekhar_retina.rds"))
cat("Saved: ", ncol(shekhar), " cells, ", nrow(shekhar), " genes\n")

cat("\n--- Aztekin Tail ---\n")
aztekin <- AztekinTailData()
saveRDS(aztekin, file.path(SC_DIR, "aztekin_tail.rds"))
cat("Saved: ", ncol(aztekin), " cells, ", nrow(aztekin), " genes\n")

# 1.14 Bunis HSPC
cat("\n--- Bunis HSPC ---\n")
bunis <- BunisHSPCData()
saveRDS(bunis, file.path(SC_DIR, "bunis_hspc.rds"))
cat("Saved: ", ncol(bunis), " cells, ", nrow(bunis), " genes\n")

cat("\n--- Chen Brain ---\n")
chen <- ChenBrainData()
saveRDS(chen, file.path(SC_DIR, "chen_brain.rds"))
cat("Saved: ", ncol(chen), " cells, ", nrow(chen), " genes\n")

cat("\n--- Darmanis Brain ---\n")
darmanis <- DarmanisBrainData()
saveRDS(darmanis, file.path(SC_DIR, "darmanis_brain.rds"))
cat("Saved: ", ncol(darmanis), " cells, ", nrow(darmanis), " genes\n")

cat("\n--- Ernst Spleen ---\n")
ernst <- ErnstSpleenData()
saveRDS(ernst, file.path(SC_DIR, "ernst_spleen.rds"))
cat("Saved: ", ncol(ernst), " cells, ", nrow(ernst), " genes\n")

cat("\n--- Fletcher Olfactory ---\n")
fletcher <- FletcherOlfactoryData()
saveRDS(fletcher, file.path(SC_DIR, "fletcher_olfactory.rds"))
cat("Saved: ", ncol(fletcher), " cells, ", nrow(fletcher), " genes\n")

# 1.19 Giladi HSPC
cat("\n--- Giladi HSPC ---\n")
giladi <- GiladiHSPCData()
saveRDS(giladi, file.path(SC_DIR, "giladi_hspc.rds"))
cat("Saved: ", ncol(giladi), " cells, ", nrow(giladi), " genes\n")

cat("\n--- He Organ Atlas ---\n")
he <- HeOrganAtlasData()
saveRDS(he, file.path(SC_DIR, "he_organ_atlas.rds"))
cat("Saved: ", ncol(he), " cells, ", nrow(he), " genes\n")

cat("\n--- Hu Cortex ---\n")
hu <- HuCortexData()
saveRDS(hu, file.path(SC_DIR, "hu_cortex.rds"))
cat("Saved: ", ncol(hu), " cells, ", nrow(hu), " genes\n")

cat("\n--- Jessa Brain Tumor ---\n")
jessa <- JessaBrainData()
saveRDS(jessa, file.path(SC_DIR, "jessa_brain_tumor.rds"))
cat("Saved: ", ncol(jessa), " cells, ", nrow(jessa), " genes\n")

# 1.23 Kolod ESC
cat("\n--- Kolodziejczyk ESC ---\n")
kolod <- KolsoeESCData()
saveRDS(kolod, file.path(SC_DIR, "kolodziejczyk_esc.rds"))
cat("Saved: ", ncol(kolod), " cells, ", nrow(kolod), " genes\n")

# 1.24 Leng ESC
cat("\n--- Leng ESC ---\n")
leng <- LengESCData()
saveRDS(leng, file.path(SC_DIR, "leng_esc.rds"))
cat("Saved: ", ncol(leng), " cells, ", nrow(leng), " genes\n")

cat("\n--- Marques Oligodendrocyte ---\n")
marques <- MarquesBrainData()
saveRDS(marques, file.path(SC_DIR, "marques_oligo.rds"))
cat("Saved: ", ncol(marques), " cells, ", nrow(marques), " genes\n")

# 1.26 Messmer ESC
cat("\n--- Messmer ESC ---\n")
messmer <- MessmerESCData()
saveRDS(messmer, file.path(SC_DIR, "messmer_esc.rds"))
cat("Saved: ", ncol(messmer), " cells, ", nrow(messmer), " genes\n")

# 1.27 Mair PBMC CITE-seq
cat("\n--- Mair PBMC ---\n")
mair <- MairPBMCData()
saveRDS(mair, file.path(SC_DIR, "mair_pbmc_cite.rds"))
cat("Saved: ", ncol(mair), " cells, ", nrow(mair), " genes\n")

# 1.28 Nestorowa HSPC
cat("\n--- Nestorowa HSPC ---\n")
nestorowa <- NestorowaHSPCData()
saveRDS(nestorowa, file.path(SC_DIR, "nestorowa_hspc.rds"))
cat("Saved: ", ncol(nestorowa), " cells, ", nrow(nestorowa), " genes\n")

cat("\n--- Pollen Brain ---\n")
pollen <- PollenGliaData()
saveRDS(pollen, file.path(SC_DIR, "pollen_brain.rds"))
cat("Saved: ", ncol(pollen), " cells, ", nrow(pollen), " genes\n")

cat("\n--- Richard T Cell ---\n")
richard <- RichardTCellData()
saveRDS(richard, file.path(SC_DIR, "richard_tcell.rds"))
cat("Saved: ", ncol(richard), " cells, ", nrow(richard), " genes\n")

# 1.31 Stoeckius PBMC Hashing
cat("\n--- Stoeckius PBMC ---\n")
stoeckius <- StoeckiusCITEseqData()
saveRDS(stoeckius, file.path(SC_DIR, "stoeckius_pbmc_hash.rds"))
cat("Saved: ", ncol(stoeckius), " cells, ", nrow(stoeckius), " genes\n")

cat("\n--- Usoskin Neuron ---\n")
usoskin <- UsoskinBrainData()
saveRDS(usoskin, file.path(SC_DIR, "usoskin_neuron.rds"))
cat("Saved: ", ncol(usoskin), " cells, ", nrow(usoskin), " genes\n")

cat("\n--- Wu Kidney ---\n")
wu <- WuKidneyData()
saveRDS(wu, file.path(SC_DIR, "wu_kidney.rds"))
cat("Saved: ", ncol(wu), " cells, ", nrow(wu), " genes\n")

cat("\n--- Xin Pancreas ---\n")
xin <- XinPancreasData()
saveRDS(xin, file.path(SC_DIR, "xin_pancreas.rds"))
cat("Saved: ", ncol(xin), " cells, ", nrow(xin), " genes\n")

cat("\n--- Zhong Prefrontal ---\n")
zhong <- ZhongPrefrontalData()
saveRDS(zhong, file.path(SC_DIR, "zhong_prefrontal.rds"))
cat("Saved: ", ncol(zhong), " cells, ", nrow(zhong), " genes\n")

cat("\n--- Zilionis Lung ---\n")
zilionis <- ZilionisLungData()
saveRDS(zilionis, file.path(SC_DIR, "zilionis_lung.rds"))
cat("Saved: ", ncol(zilionis), " cells, ", nrow(zilionis), " genes\n")

# =============================================================================
# =============================================================================

cat("\n========== : TENxPBMCData ==========\n")

# 2.1 PBMC 3k
cat("\n--- PBMC 3k ---\n")
pbmc3k <- TENxPBMCData("pbmc3k")
saveRDS(pbmc3k, file.path(SC_DIR, "tenx_pbmc3k.rds"))
cat("Saved: ", ncol(pbmc3k), " cells\n")

# 2.2 PBMC 4k
cat("\n--- PBMC 4k ---\n")
pbmc4k <- TENxPBMCData("pbmc4k")
saveRDS(pbmc4k, file.path(SC_DIR, "tenx_pbmc4k.rds"))
cat("Saved: ", ncol(pbmc4k), " cells\n")

# 2.3 PBMC 8k
cat("\n--- PBMC 8k ---\n")
pbmc8k <- TENxPBMCData("pbmc8k")
saveRDS(pbmc8k, file.path(SC_DIR, "tenx_pbmc8k.rds"))
cat("Saved: ", ncol(pbmc8k), " cells\n")

# 2.4 PBMC 6k
cat("\n--- PBMC 6k ---\n")
pbmc6k <- TENxPBMCData("pbmc6k")
saveRDS(pbmc6k, file.path(SC_DIR, "tenx_pbmc6k.rds"))
cat("Saved: ", ncol(pbmc6k), " cells\n")

# 2.5 Frozen PBMC Donor A
cat("\n--- Frozen PBMC Donor A ---\n")
frozen_a <- TENxPBMCData("frozen_pbmc_donor_a")
saveRDS(frozen_a, file.path(SC_DIR, "tenx_frozen_a.rds"))
cat("Saved: ", ncol(frozen_a), " cells\n")

# 2.6 Frozen PBMC Donor B
cat("\n--- Frozen PBMC Donor B ---\n")
frozen_b <- TENxPBMCData("frozen_pbmc_donor_b")
saveRDS(frozen_b, file.path(SC_DIR, "tenx_frozen_b.rds"))
cat("Saved: ", ncol(frozen_b), " cells\n")

# 2.7 Frozen PBMC Donor C
cat("\n--- Frozen PBMC Donor C ---\n")
frozen_c <- TENxPBMCData("frozen_pbmc_donor_c")
saveRDS(frozen_c, file.path(SC_DIR, "tenx_frozen_c.rds"))
cat("Saved: ", ncol(frozen_c), " cells\n")

# =============================================================================
# =============================================================================

cat("\n========== : muscData ==========\n")

cat("\n--- Kang PBMC 8vs8 ---\n")
kang <- Kang18_8vs8()
saveRDS(kang, file.path(SC_DIR, "kang_pbmc_8vs8.rds"))
cat("Saved: ", ncol(kang), " cells, ", length(unique(kang$ind)), " donors\n")

# =============================================================================
# =============================================================================

cat("\n========== : ExperimentHub ==========\n")

eh <- ExperimentHub()

# 4.1 Tabula Muris
cat("\n--- Tabula Muris (FACS) ---\n")
tm_query <- query(eh, "TabulaMuris")
tm_facs <- eh[["EH1617"]]  # Tabula Muris FACS
saveRDS(tm_facs, file.path(SC_DIR, "tabula_muris_facs.rds"))
cat("Saved Tabula Muris FACS\n")

# =============================================================================
# =============================================================================

cat("\n========== Bulk : GTEx (recount3) ==========\n")

human_projects <- available_projects()
gtex_projects <- subset(human_projects, file_source == "gtex")

gtex_tissues <- unique(gtex_projects$project)
cat("GTEx tissues: ", length(gtex_tissues), "\n")

sapply(gtex_tissues, function(tissue) {
    cat("\n--- GTEx:", tissue, "---\n")
    
    result <- tryCatch({
        info <- subset(gtex_projects, project == tissue)[1, ]
        rse <- create_rse(info)
        
        counts <- assay(rse, "raw_counts")
        metadata <- as.data.frame(colData(rse))
        
        saveRDS(list(counts = counts, metadata = metadata), 
                file.path(BULK_DIR, paste0("gtex_", tolower(tissue), ".rds")))
        cat("Saved: ", ncol(counts), " samples, ", nrow(counts), " genes\n")
        "OK"
    }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        "FAILED"
    })
    NULL
})

# =============================================================================
# =============================================================================

cat("\n========== Bulk : TCGA ==========\n")

tcga_cancers <- c(
    "ACC",   # Adrenocortical carcinoma
    "BLCA",  # Bladder Urothelial Carcinoma
    "BRCA",  # Breast invasive carcinoma
    "CESC",  # Cervical squamous cell carcinoma
    "CHOL",  # Cholangiocarcinoma
    "COAD",  # Colon adenocarcinoma
    "DLBC",  # Diffuse Large B-cell Lymphoma
    "ESCA",  # Esophageal carcinoma
    "GBM",   # Glioblastoma multiforme
    "HNSC",  # Head and Neck squamous cell carcinoma
    "KICH",  # Kidney Chromophobe
    "KIRC",  # Kidney renal clear cell carcinoma
    "KIRP",  # Kidney renal papillary cell carcinoma
    "LAML",  # Acute Myeloid Leukemia
    "LGG",   # Brain Lower Grade Glioma
    "LIHC",  # Liver hepatocellular carcinoma
    "LUAD",  # Lung adenocarcinoma
    "LUSC",  # Lung squamous cell carcinoma
    "MESO",  # Mesothelioma
    "OV",    # Ovarian serous cystadenocarcinoma
    "PAAD",  # Pancreatic adenocarcinoma
    "PCPG",  # Pheochromocytoma and Paraganglioma
    "PRAD",  # Prostate adenocarcinoma
    "READ",  # Rectum adenocarcinoma
    "SARC",  # Sarcoma
    "SKCM",  # Skin Cutaneous Melanoma
    "STAD",  # Stomach adenocarcinoma
    "TGCT",  # Testicular Germ Cell Tumors
    "THCA",  # Thyroid carcinoma
    "THYM",  # Thymoma
    "UCEC",  # Uterine Corpus Endometrial Carcinoma
    "UCS",   # Uterine Carcinosarcoma
    "UVM"    # Uveal Melanoma
)

sapply(tcga_cancers, function(cancer) {
    cat("\n--- TCGA:", cancer, "---\n")
    
    result <- tryCatch({
        mae <- curatedTCGAData(cancer, "RNASeq2Gene", dry.run = FALSE, version = "2.0.1")
        
        expr <- assay(mae[[1]])
        
        sample_codes <- TCGAbarcode(colnames(expr), sample = TRUE)
        sample_types <- substr(sample_codes, 1, 2)
        
        tumor_n <- sum(sample_types == "01")
        normal_n <- sum(sample_types == "11")
        
        saveRDS(list(
            expression = expr,
            sample_types = sample_types,
            tumor_n = tumor_n,
            normal_n = normal_n
        ), file.path(BULK_DIR, paste0("tcga_", tolower(cancer), ".rds")))
        
        cat("Saved: Tumor=", tumor_n, ", Normal=", normal_n, "\n")
        "OK"
    }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        "FAILED"
    })
    NULL
})

# =============================================================================
# =============================================================================

cat("\n========== Bulk : GEO ==========\n")

geo_datasets <- c(
    "GSE147507",
    "GSE68086",
    "GSE48350",
    "GSE22260",
    "GSE62944",   # TCGA Pan-Cancer
    "GSE81861",
    "GSE84133",
    "GSE63818",
    "GSE75140",
    "GSE36552"
)

sapply(geo_datasets, function(gse) {
    cat("\n--- GEO:", gse, "---\n")
    
    result <- tryCatch({
        gset <- getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE)
        
        expr <- exprs(gset[[1]])
        pdata <- pData(gset[[1]])
        
        saveRDS(list(expression = expr, phenotype = pdata), 
                file.path(BULK_DIR, paste0("geo_", tolower(gse), ".rds")))
        
        cat("Saved: ", ncol(expr), " samples, ", nrow(expr), " features\n")
        "OK"
    }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        "FAILED"
    })
    NULL
})

# =============================================================================
# =============================================================================

cat("\n========== ==========\n")

sc_files <- list.files(SC_DIR, pattern = "\\.rds$", full.names = FALSE)
bulk_files <- list.files(BULK_DIR, pattern = "\\.rds$", full.names = FALSE)

cat("\n:", length(sc_files), "\n")
cat(paste(" -", sc_files, collapse = "\n"), "\n")

cat("\nBulk :", length(bulk_files), "\n")
cat(paste(" -", bulk_files, collapse = "\n"), "\n")

cat("\n:", DATA_DIR, "\n")
cat("\n========== ==========\n")
