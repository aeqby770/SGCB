PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
STORY_ROOT <- normalizePath(paste0(PROJECT_ROOT, "/biological_story_workspace"), winslash = "/", mustWork = TRUE)
STORY_OUTPUT_ROOT <- file.path(STORY_ROOT, "story_outputs")
dir.create(STORY_OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

STORY_DEFAULT_CASES <- c("tcga_kirc", "pb_kang18", "proteomics_cptac_ccrcc")
STORY_FDR_CUTOFF <- 0.05
STORY_TOP_FEATURES_PER_CHANNEL <- 6L
STORY_TOP_ENRICHMENT_TERMS <- 10L
STORY_SIGNATURE_TOP_K <- 30L

# Group colors for visualization (must cover group_levels across all cases)
BIO_GROUP_COLORS <- c(
  normal     = "#4C72B0",  tumor      = "#C44E52",
  ctrl       = "#4C72B0",  stim       = "#C44E52",
  control    = "#4C72B0",  treated    = "#C44E52",
  low_grade  = "#4C72B0",  high_grade = "#C44E52"
)

STORY_CASES <- list(
  tcga_kirc = list(
    case_id = "tcga_kirc",
    title = "TCGA-KIRC tumor vs normal",
    modality = "bulk",
    species = "human",
    input_path = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/04_tcga/tcga_kirc_tumor_vs_normal.rds"),
    cached_result_path = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/sgcb/results/tcga__kirc_full.rds"),
    output_dir = file.path(STORY_OUTPUT_ROOT, "tcga_kirc"),
    group_levels = c("normal", "tumor"),
    use_cached_default = TRUE,
    representative_prefix = c(
      "DE-only" = "de",
      "DV-only" = "dv",
      "DG-only" = "dg",
      "overlap" = "ov"
    )
  ),
  pb_kang18 = list(
    case_id = "pb_kang18",
    title = "PB-kang18 pseudobulk stimulation",
    modality = "pseudobulk",
    species = "human",
    input_path = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/data/bulk/05_pseudobulk/kang18_pseudobulk.rds"),
    cached_result_path = paste0(PROJECT_ROOT, "/benchmark_v2/01_bulk_rnaseq/sgcb/results/pseudobulk__kang18_full.rds"),
    output_dir = file.path(STORY_OUTPUT_ROOT, "pb_kang18"),
    group_levels = NULL,
    use_cached_default = TRUE,
    representative_prefix = c(
      "DE-only" = "de",
      "DV-only" = "dv",
      "DG-only" = "dg",
      "overlap" = "ov"
    )
  ),
  proteomics_cptac_ccrcc = list(
    case_id = "proteomics_cptac_ccrcc",
    title = "CPTAC-CCRCC proteome low grade vs high grade",
    modality = "proteomics",
    species = "human",
    input_path = paste0(PROJECT_ROOT, "/biological_story_workspace/proteomics/CPTAC_CCRCC/PDC000127_proteome_log2ratio.csv"),
    phospho_path = paste0(PROJECT_ROOT, "/biological_story_workspace/proteomics/CPTAC_CCRCC/PDC000128_phosphoproteome_log2ratio.csv"),
    clinical_path = paste0(PROJECT_ROOT, "/biological_story_workspace/proteomics/CPTAC_CCRCC/PDC000127_clinical_metadata.csv"),
    cached_result_path = NULL,
    output_dir = file.path(STORY_OUTPUT_ROOT, "proteomics_cptac_ccrcc"),
    group_levels = c("low_grade", "high_grade"),
    grade_low = c("G1", "G2"),
    grade_high = c("G3", "G4"),
    use_cached_default = FALSE,
    representative_prefix = c(
      "DE-only" = "de",
      "DV-only" = "dv",
      "DG-only" = "dg",
      "overlap" = "ov"
    )
  )
)

story_case_ids <- function() {
  names(STORY_CASES)
}

story_get_case <- function(case_id) {
  STORY_CASES[[case_id]]
}
