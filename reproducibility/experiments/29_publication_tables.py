from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import beta

CONFIG = {"paths": {"root": Path(sys.argv[1]),
    "comparator_results":Path(sys.argv[1])/"tmp/revision_0.5.5_integrated/results"}, "threads": int(os.environ["SGCB_THREADS"]), "plot": {}, "cohorts": ["sle_classical_monocytes_2021","masld_liver_2023"],
    "cohort_labels": {"sle_classical_monocytes_2021":"SLE monocytes","masld_liver_2023":"MASLD liver"},
    "mean_fit_version": "0.5.7",
    "publication_metadata_columns": ['run_id', 'fit_started', 'native_source_id', 'mean_inference_source_id', 'input_id', 'version', 'release_version', 'fit_package_version', 'fit_dataset', 'fit_version', 'sgcb_fit_version', 'sgcb_reuse_reason', 'reuse_reason', 'result_label'],
    "publication_label_replacements": {'^sgcb_\\d+\\.\\d+\\.\\d+_': '', '^current_kirc$': 'KIRC', '_streamed$': '', '^current arithmetic-mean public API$': 'SGCB arithmetic-mean inference', '^unchanged archived matched-input fit$': 'Matched-input comparator fit', '^CLR of archived counts plus 0\\.5$': 'CLR of counts plus 0.5', '^current public API DV component$': 'API DV component', '^current Gamma DV refit$': 'Gamma-restricted DV fit', '^current lognormal DV refit$': 'Lognormal-restricted DV fit', '^complete retained output$': 'available for all partitions', '^partial retained output$': 'available for some partitions', '^archived comparator context; no current population fit$': 'Gene-wise GG and comparator analysis in GTEx brain tissues', '^historical gene-wise covariance diagnostic$': 'gene-wise covariance diagnostic'}}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
plt.rcParams.update(CONFIG["plot"])
ROOT = CONFIG["paths"]["root"]
VERSION = (ROOT / "software" / "SGCB" / "DESCRIPTION").read_text(encoding="utf-8").split("Version: ",1)[1].splitlines()[0]
R = ROOT / "results"
OUT = ROOT / "manuscript"

# Population-sharing results retain every fixed design and truth domain.
population_calibration = pd.read_csv(R / "publication_population_family_summary.csv")
population_calibration = population_calibration.loc[population_calibration.family.eq("GG")].copy()
population_scenarios = pd.read_csv(R / "publication_population_family_scenarios.csv")
population_scenarios = population_scenarios.loc[population_scenarios.family.eq("GG")].copy()
population_calibration["calls"] = (population_calibration.TP + population_calibration.FP +
    population_calibration.calls_with_undefined_truth)
current_api_metrics = pd.read_csv(R/"publication_current_api_metrics.csv")
current_api_scenarios = pd.read_csv(R/"publication_current_api_scenarios.csv")
current_api_metrics["version"] = VERSION
current_api_scenarios["version"] = VERSION

population_calibration.to_csv(R/"publication_population_calibration.csv",index=False,encoding="utf-8")
population_scenarios.to_csv(R/"publication_population_scenarios.csv",index=False,encoding="utf-8")

population_family_summary = pd.read_csv(R / "publication_population_family_summary.csv")
population_family_scenarios = pd.read_csv(R / "publication_population_family_scenarios.csv")
population_family_operating = pd.read_csv(R / "publication_population_family_operating_curves.csv")
# Independent applications use the same primary DE/DV definitions.
COHORTS = CONFIG["cohorts"]
bio_prefix = str(R)+f"/sgcb_{VERSION}_"+pd.Series(COHORTS)
population_biology_genes = pd.concat(map(pd.read_csv,bio_prefix+"_biological_results.csv"),
    keys=COHORTS,names=["dataset","record"]).reset_index().drop(columns="record")
population_biology_genes["conventional_mean_available"] = population_biology_genes.any_conventional_DE.notna()
population_biology_genes["primary_DV_only"] = (
    population_biology_genes.population_DV & ~population_biology_genes.DE)
population_biology_genes["primary_adjusted_DV_only"] = (
    population_biology_genes.primary_DV_only & population_biology_genes.population_adjusted_retained)
population_biology_genes["any_mean_negative"] = (
    population_biology_genes.population_DV & ~population_biology_genes.DE &
    ~population_biology_genes.any_conventional_DE.astype("boolean").fillna(False)).where(
        population_biology_genes.conventional_mean_available)
population_biology_genes["any_mean_negative_retained"] = (
    population_biology_genes.population_DV & ~population_biology_genes.DE &
    ~population_biology_genes.any_conventional_DE.astype("boolean").fillna(False) &
    population_biology_genes.population_adjusted_retained).where(
        population_biology_genes.conventional_mean_available)
population_biology = population_biology_genes.groupby("dataset",sort=False).agg(
    n_genes=("gene_id","size"),SGCB_DE=("DE","sum"),SGCB_DV=("population_DV","sum"),
    DV_available=("population_available","sum"),DV_only=("primary_DV_only","sum"),
    retained_DV_only=("primary_adjusted_DV_only","sum"),
    conventional_mean_available=("conventional_mean_available","all")).reset_index()
mean_sensitivity = population_biology_genes.groupby("dataset",sort=False)[
    ["any_mean_negative","any_mean_negative_retained"]].sum(min_count=1).reset_index()
population_biology = population_biology.merge(mean_sensitivity,on="dataset",validate="one_to_one")
bio_sizes = pd.concat(map(pd.read_csv,bio_prefix+"_summary.csv"),ignore_index=True)[
    ["dataset","n_ctrl","n_treat"]].rename(columns={"n_ctrl":"n_control","n_treat":"n_case"})
bio_sizes["dataset"] = bio_sizes.dataset.str.removeprefix(f"sgcb_{VERSION}_")
population_biology = population_biology.merge(bio_sizes,on="dataset",validate="one_to_one")
sle_covariates = pd.read_csv(R/f"sgcb_{VERSION}_sle_classical_monocytes_2021_covariate_summary.csv").rename(
    columns={"n_genes":"features","calls":"moment_calls","SGCB_selected":"population_selected",
        "selected_sign_retained":"selected_direction_retained"})
population_covariates = pd.concat([
    sle_covariates.assign(dataset=COHORTS[0]),
    pd.read_csv(R/f"sgcb_{VERSION}_masld_liver_2023_covariate_summary.csv").assign(dataset=COHORTS[1])],
    ignore_index=True)
population_biology.to_csv(R/"publication_population_biology.csv",index=False,encoding="utf-8")
population_covariates.to_csv(R/"publication_population_covariates.csv",index=False,encoding="utf-8")

null = pd.concat(map(pd.read_csv,
    CONFIG["paths"]["comparator_results"].glob("readcount_independent_null_comparator_*_metrics.csv")),ignore_index=True).rename(
        columns={"genes":"n_common_genes","available":"n_valid_p","has_call":"any_discovery"})
null = null.loc[~null.method.isin(["SGCB","Gamma","lognormal"])].copy()
null["method"] = null.method.replace({"diffVar":"DiffVar","Brown_Forsythe":"Brown-Forsythe"})
null_sample_sizes = null[["split","n_A","n_B"]].drop_duplicates()
current_null = pd.read_csv(R/f"sgcb_{VERSION}_readcount_null_metrics.csv")
null = pd.concat([null,current_null],ignore_index=True).sort_values(["method","target","split"])
null["positive_fraction"] = null.calls/null.n_common_genes
null_summary = null.groupby(["method","target"],sort=False).agg(
    partitions=("split","nunique"),with_discovery=("any_discovery","sum"),
    mean_calls=("calls","mean"),min_calls=("calls","min"),max_calls=("calls","max"),
    min_common_genes=("n_common_genes","min"),max_common_genes=("n_common_genes","max"),
    mean_raw_p005=("raw_p005_full","mean"),mean_conditional_coverage=("null_coverage","mean")).reset_index()
null_summary["partition_fraction_with_calls"] = null_summary.with_discovery/null_summary.partitions
external_methods = pd.DataFrame({
    "method":["DESeq2","edgeR","limma","DiffVar","clrDV","GAMLSS","MDSeq","Brown-Forsythe"],
    "target":["Mean expression","Mean expression","Mean expression","Log spread",
        "CLR dispersion","Log spread","Count dispersion","Log spread"]})
null_availability = pd.concat([current_null[["method","target"]].drop_duplicates(),
    external_methods],ignore_index=True).merge(
        null_summary[["method","target","partitions"]],on=["method","target"],how="left",
        validate="one_to_one").rename(columns={"partitions":"completed_partitions"})
null_availability["completed_partitions"] = null_availability.completed_partitions.fillna(0).astype(int)
null_availability["fixed_partitions"] = 10
null_availability["missing_partitions"] = 10-null_availability.completed_partitions
null_availability["status"] = np.select([
    null_availability.completed_partitions.eq(10),null_availability.completed_partitions.gt(0)],
    ["available for all partitions","available for some partitions"],default="matched output unavailable")


# Gene-wise GG and native comparator summaries for GTEx brain tissues.
brain = pd.read_csv(R / "current_brain_common_summary.csv")
brain["method"] = brain.method.replace({"SGCB":"Gene-wise GG"})
brain_authors = pd.concat([
    pd.read_csv(R / "brain_cerebellum_vs_cortex_diffVar_summary.csv"),
    pd.read_csv(R / "brain_cerehemi_vs_cortex_diffVar_summary.csv"),
    pd.read_csv(R / "brain_hippo_vs_hypo_diffVar_summary.csv"),
    pd.read_csv(R / "brain_cerebellum_vs_cortex_clrDV_summary.csv"),
    pd.read_csv(R / "brain_cerehemi_vs_cortex_clrDV_summary.csv"),
    pd.read_csv(R / "brain_hippo_vs_hypo_clrDV_summary.csv")],ignore_index=True)
brain_authors["method"] = brain_authors.method.replace({"diffVar":"DiffVar"})
brain = pd.concat([brain,
    brain_authors.assign(contrast=brain_authors.dataset.str.removeprefix("brain_"),
        calls=brain_authors.BH_calls,n_common_genes=brain_authors.genes),
    brain_authors.loc[brain_authors.method.eq("clrDV")].assign(method="clrDV (author BY)",
        contrast=brain_authors.dataset.str.removeprefix("brain_"),
        calls=brain_authors.author_calls,n_common_genes=brain_authors.genes)],ignore_index=True)
brain = brain.rename(columns={"overlap_SGCB":"overlap_gene_wise_GG",
    "jaccard_SGCB":"jaccard_gene_wise_GG","reference_SGCB_calls":"reference_gene_wise_GG_calls"})
brain["analysis_role"] = "Gene-wise GG and comparator analysis in GTEx brain tissues"
qpcr = pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_seqc_qpcr_metrics.csv")
qpcr["sgcb_fit_version"] = pd.Series(CONFIG["mean_fit_version"], index=qpcr.index).where(qpcr.method.eq("SGCB"))
qpcr["sgcb_used_endpoint"] = pd.Series("arithmetic_mean_DE", index=qpcr.index).where(qpcr.method.eq("SGCB"))
qpcr["sgcb_reuse_reason"] = pd.Series("saved mean-endpoint fit on the same inputs", index=qpcr.index).where(qpcr.method.eq("SGCB"))
sim = pd.read_csv(R / "current_simulation_metrics.csv").assign(method="Gene-wise GG",inference_role="diagnostic log-mean / CV2 baseline")
sim_scenario = pd.read_csv(R / "current_simulation_scenarios.csv").assign(method="Gene-wise GG",inference_role="diagnostic log-mean / CV2 baseline")
proteins_legacy = pd.concat([
    pd.read_csv(R / "protein_ecoli_lfq_common_metrics.csv").assign(dataset="E. coli LFQ"),
    pd.read_csv(R / "protein_sgsds_ratio2_common_metrics.csv").assign(dataset="SGSDS ratio 2"),
    pd.read_csv(R / "protein_sgsds_ratio2.5_common_metrics.csv").assign(dataset="SGSDS ratio 2.5"),
    pd.read_csv(R / "protein_cptac_ups1_common_metrics.csv").assign(dataset="CPTAC UPS1")], ignore_index=True)
protein_matched = pd.concat([
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_protein_ecoli_lfq_metrics.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_protein_sgsds_ratio2_metrics.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_protein_sgsds_ratio2.5_metrics.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_protein_cptac_ups1_metrics.csv")], ignore_index=True)
protein_matched["sgcb_fit_version"] = pd.Series(CONFIG["mean_fit_version"], index=protein_matched.index).where(protein_matched.method.eq("SGCB"))
protein_matched["sgcb_used_endpoint"] = pd.Series("arithmetic_mean_DE", index=protein_matched.index).where(protein_matched.method.eq("SGCB"))
protein_matched["score_source"] = protein_matched.score_source.replace({"current arithmetic-mean public API":"SGCB arithmetic-mean inference", "unchanged archived matched-input fit":"Matched-input comparator fit"})
protein_matched["F1"] = 2 * protein_matched.TP / (
    protein_matched.spike_proteins + protein_matched.TP + protein_matched.FP)
proteins = protein_matched[["dataset", "method", "normalization", "features",
    "spike_proteins", "null_proteins", "TP", "FP", "power", "FDP", "F1",
    "raw_null_rejection", "score_source", "sgcb_fit_version", "sgcb_used_endpoint"]].copy()
proteins.to_csv(R / "publication_proteomics_matched.csv", index=False, encoding="utf-8")

nb = pd.concat([
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_nb_sim_n3_pDE10_ef1.5_summary.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_nb_sim_n3_pDE10_ef2.0_summary.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_nb_sim_n5_pDE10_ef1.5_summary.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_nb_sim_n5_pDE10_ef2.0_summary.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_nb_sim_n10_pDE10_ef1.5_summary.csv"),
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_nb_sim_n10_pDE10_ef2.0_summary.csv")], ignore_index=True)
nb["fit_version"] = CONFIG["mean_fit_version"]
nb["used_endpoint"] = "arithmetic_mean_DE"
cv_adjustment = pd.read_csv(R / "publication_CV2_adjustment.csv")
cv_kirc = pd.read_csv(R / f"sgcb_{VERSION}_kirc_covariate_genes.csv")
cv_kirc = cv_kirc.assign(adjusted_significant=cv_kirc.adjusted_qvalue.lt(.05),
    direction_retained=np.sign(cv_kirc.adjusted_log2_cv2_ratio).eq(np.sign(cv_kirc.primary_GG_effect)))
cv_by_channel = cv_kirc.groupby("channel",sort=False).agg(
    genes=("gene_id","size"),adjusted_available=("adjusted_pvalue","count"),
    adjusted_significant=("adjusted_significant","sum"),direction_retained=("direction_retained","sum")).reset_index()
covariance_diagnosis = pd.concat([
    pd.read_csv(R / "penalized_covariance_abundance_n50_50_summary.csv").assign(design="50:50",input="Abundance"),
    pd.read_csv(R / "penalized_covariance_abundance_n10_50_summary.csv").assign(design="10:50",input="Abundance"),
    pd.read_csv(R / "penalized_covariance_abundance_n10_10_summary.csv").assign(design="10:10",input="Abundance"),
    pd.read_csv(R / "penalized_covariance_counts_n50_50_summary.csv").assign(design="50:50",input="Liver null 1")],ignore_index=True)
survival = pd.concat([pd.read_csv(R / f"sgcb_{VERSION}_kirc_heldout_survival.csv").assign(cohort="TCGA held-out"),
                      pd.read_csv(R / f"sgcb_{VERSION}_kirc_external_survival.csv").assign(cohort="E-MTAB-1980")], ignore_index=True)
alternative = pd.concat([pd.read_csv(R / f"sgcb_{VERSION}_kirc_alternative_scores_survival.csv").assign(cohort="TCGA held-out"),
                         pd.read_csv(R / f"sgcb_{VERSION}_kirc_external_alternative_scores_survival.csv").assign(cohort="E-MTAB-1980")], ignore_index=True)
distribution_scores = pd.concat([
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_variability_heldout_survival.csv"),
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_variability_external_survival.csv")], ignore_index=True)
ph = pd.concat([pd.read_csv(R / f"sgcb_{VERSION}_kirc_heldout_proportional_hazards.csv").assign(cohort="TCGA held-out"),
                pd.read_csv(R / f"sgcb_{VERSION}_kirc_external_proportional_hazards.csv").assign(cohort="E-MTAB-1980"),
                pd.read_csv(R / f"sgcb_{VERSION}_kirc_alternative_scores_PH.csv").rename(columns={"score":"signature"}).assign(cohort="TCGA held-out"),
                pd.read_csv(R / f"sgcb_{VERSION}_kirc_external_alternative_scores_PH.csv").rename(columns={"score":"signature"}).assign(cohort="E-MTAB-1980")], ignore_index=True)
study = pd.DataFrame({
    "Dataset": ["GG / misspecification simulation","Negative-binomial simulations","GTEx liver null","GTEx brain: cerebellum / cortex","GTEx brain: hemisphere / cortex","GTEx brain: hippocampus / hypothalamus","SEQC reference RNA","TCGA-KIRC discovery","TCGA-KIRC held-out survival","E-MTAB-1980 survival","Kang pseudobulk","Protein spike-ins"],
    "Group sizes / units": ["50:50; 10:10; 10:50","3:3, 5:5, 10:10","50:50 per partition","50:50; 93 donors","50:50; 92 donors","50:50; 84 donors","5:5 technical preparations","59 normals; 80 tumors; 139 patients","392 patients; 131 deaths","101 patients; 23 deaths","8 paired donors; 16 profiles","3:3 per contrast"],
    "Features / contrasts": ["5,400; nine scenarios","10,000; six contrasts","49,061–50,443 common genes; ten partitions","52,454 common genes","52,306 common genes","52,094 common genes","18,615 SGCB-tested; 683 qPCR-common","20,031 genes","385 clinically complete","99 clinically complete","10,385 genes","Four contrasts; 704–4,508 complete proteins"]})
study.loc[2,"Group sizes / units"] = (
    f"{null_sample_sizes.n_A.min():d}–{null_sample_sizes.n_A.max():d} / "
    f"{null_sample_sizes.n_B.min():d}–{null_sample_sizes.n_B.max():d} independent donors")
study.loc[2,"Features / contrasts"] = (
    f"{null.n_common_genes.min():,}–{null.n_common_genes.max():,} common genes; ten liver partitions")
study = pd.concat([study,pd.DataFrame({
    "Dataset":["Observed-zero calibration","SLE classical monocytes (2021)","Metabolic liver disease (2023)"],
    "Group sizes / units":["50:50; 10:10; 10:50","24 controls; 64 patients; first visits","69 controls; 299 patients"],
    "Features / contrasts":["4,000; 20 scenarios; 3,800 finite-CV2 truths",
        "23,972 genes; integer read counts","19,787 genes; deposited FPKM"]})],ignore_index=True)
new_biology = population_biology
new_biology_genes = population_biology_genes
new_covariates = population_covariates
new_enrichment = pd.concat([
    pd.read_csv(R/f"sgcb_{VERSION}_sle_classical_monocytes_2021_enrichment.csv").assign(dataset=COHORTS[0]),
    pd.read_csv(R/f"sgcb_{VERSION}_masld_liver_2023_enrichment.csv").assign(dataset=COHORTS[1])],ignore_index=True)
table4 = new_biology.assign(Cohort=new_biology.dataset.map(CONFIG["cohort_labels"]),
    size=new_biology.n_control.astype(str)+" / "+new_biology.n_case.astype(str))[
        ["Cohort","size","n_genes","SGCB_DV","DV_only","retained_DV_only"]]
table4.columns = ["Cohort","Control / case","Genes","DV","DV-only","Adjusted DV-only"]
new_biology.to_csv(R/"publication_new_biology.csv",index=False,encoding="utf-8")
table2 = qpcr.loc[pd.to_numeric(qpcr.qpcr_abs_log2FC_threshold,errors="coerce").eq(1) &
    qpcr.reference_estimand.eq("geometric_mean")].sort_values("AUC", ascending=False).copy()
table2["AUC (95% CI)"] = table2.AUC.round(3).astype(str) + " (" + table2.AUC_low95.round(3).astype(str) + "–" + table2.AUC_high95.round(3).astype(str) + ")"
table2 = table2[["method","AUC (95% CI)","average_precision","F1","precision","recall"]]
table2.columns = ["Method","AUC (95% CI)","AP","F1","Precision","Recall"]
table3 = pd.DataFrame({
    "Method":null_summary.method,"Target":null_summary.target,
    "With calls":null_summary.with_discovery.astype(int).astype(str)+"/"+null_summary.partitions.astype(int).astype(str),
    "Calls: mean (range)":null_summary.mean_calls.round(1).astype(str)+" ("+
        null_summary.min_calls.astype(int).astype(str)+"–"+null_summary.max_calls.astype(int).astype(str)+")",
    "Genes":null_summary.min_common_genes.astype(int).astype(str)+"–"+
        null_summary.max_common_genes.astype(int).astype(str)})
missing_null_methods = null_availability.loc[null_availability.completed_partitions.eq(0)]
table3 = pd.concat([table3,pd.DataFrame({"Method":missing_null_methods.method,
    "Target":missing_null_methods.target,"With calls":"unavailable",
    "Calls: mean (range)":"—","Genes":"—"})],ignore_index=True)
study = study.loc[~study.Dataset.str.startswith("GTEx brain:")].copy()
study_lines = study.to_csv(sep="|",index=False).splitlines()
study_markdown = "\n".join(study_lines[:1]+["|".join(["---"]*len(study.columns))]+study_lines[1:])
table2_lines = table2.to_csv(sep="|",index=False,float_format="%.3f").splitlines()
table2_markdown = "\n".join(table2_lines[:1]+["|".join(["---"]*len(table2.columns))]+table2_lines[1:])
table3_lines = table3.to_csv(sep="|",index=False).splitlines()
table3_markdown = "\n".join(table3_lines[:1]+["|".join(["---"]*len(table3.columns))]+table3_lines[1:])
table4_lines = table4.to_csv(sep="|",index=False).splitlines()
table4_markdown = "\n".join(table4_lines[:1]+["|".join(["---"]*len(table4.columns))]+table4_lines[1:])
table_text = "## Tables\n\n### Table 1. Evaluation datasets and biological units\n\n" + study_markdown
table_text += "\n\n### Table 2. External TaqMan evaluation on 683 common genes\n\n" + table2_markdown
table_text += "\n\nThe primary qPCR effect threshold is |log2(B/A)| ≥ 1 using the geometric-mean reference. The matched arithmetic-mean sensitivity uses the same genes and thresholds in Table S5. AP denotes average precision. AUC intervals are approximate DeLong intervals; all F1 intervals and alternative thresholds are in Table S5. Adjusted p-values retain each method's full tested family.\n\n### Table 3. Discoveries in the ten fixed GTEx liver null partitions\n\n" + table3_markdown
table_text += "\n\nAll ten liver partitions use integer read counts reconstructed at 40 million reads, removal of every repeated-donor record, and the shared expression filter and TMM normalization. The sample sizes and exact common gene families appear in Table S2. Methods retain their stated abundance-CV2 or log-spread targets. Counts report observed false discoveries and the fraction of partitions with any call; overlapping partitions are described without treating them as independent Bernoulli trials. Matched retained outputs cover all ten partitions for DESeq2, edgeR, limma, DiffVar, GAMLSS and Brown–Forsythe, and partition 1 for MDSeq. The clrDV result is unavailable for these inputs; Table S2 lists comparator coverage by partition.\n"
table_text += "\n\n### Table 4. Independent disease applications and covariate sensitivity\n\n" + table4_markdown
table_text += ("\n\nDV-only excludes SGCB DE in the same complete gene family. Adjusted DV-only counts same-direction significant conditional-moment CV2 contrasts among that set, after age and sex for SLE or clinical severity and composition proxies for liver. Conventional mean comparisons are available for liver only: "
    + str(int(new_biology.set_index("dataset").loc["masld_liver_2023","any_mean_negative"]))
    + " DV genes lack SGCB and conventional mean evidence, with "
    + str(int(new_biology.set_index("dataset").loc["masld_liver_2023","any_mean_negative_retained"]))
    + " retaining same-direction full-adjusted evidence. SLE mean comparisons are SGCB-only. Tables S24 and S26 provide gene-level results, covariate designs, sample metadata and unavailable comparisons. Each profile represents one individual.\n")
(OUT / "tables.md").write_text(table_text, encoding="utf-8")
with pd.ExcelWriter(OUT / "Supplementary_Tables.xlsx", engine="openpyxl") as writer:
    study.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S1_Study_design", index=False)
    null.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S2_Null_partitions", index=False)
    null_summary.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S2_Null_aggregate", index=False)
    null_availability.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S2_Method_availability",index=False)
    sim.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S3_Simulation", index=False)
    sim_scenario.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S3_Scenario", index=False)
    pd.read_csv(R / "publication_CV2_comparison.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S3_Same_CV2_tests", index=False)
    pd.read_csv(R / "current_CV2_estimation.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S4_Moment_estimation", index=False)
    qpcr.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S5_qPCR", index=False)
    current_api_metrics.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S27_Endpoint_metrics", index=False)
    current_api_scenarios.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S27_Scenario_metrics", index=False)
    brain.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S6_Gene_wise_brain", index=False)
    cv_adjustment.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S7_Conditional_CV2", index=False)
    cv_by_channel.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S7_KIRC_CV2_by_channel", index=False)
    proteins.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S8_Proteomics", index=False)
    proteins_legacy.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S8_Preprocess_sensitivity", index=False)
    nb.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S9_NB_simulations", index=False)
    pd.read_csv(R / f"sgcb_{CONFIG['mean_fit_version']}_pseudobulk_kang18_paired_summary.csv").assign(
        fit_version=CONFIG["mean_fit_version"],used_endpoint="arithmetic_mean_DE").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S10_Pseudobulk", index=False)
    survival.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S11_Survival", index=False)
    alternative.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S11_Score_sensitivity", index=False)
    distribution_scores.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S11_Distribution_scores", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_variability_discovery_score_parameters.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
        writer, sheet_name="S11_Distribution_parameters", index=False)
    ph.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S12_PH_diagnostics", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_external_coverage.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S13_External_coverage", index=False)
    pd.read_csv(R/"publication_current_resources.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S14_Complete_API",index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_summary.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S15_KIRC_channels", index=False)
    pd.read_csv(R / "current_kirc_summary.csv").drop(columns=["run_id", "fit_started", "native_source_id", "input_id"]).drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S15_Gene_wise_KIRC", index=False)
    pd.read_csv(R / "current_zero_fraction_fit_summary.csv",keep_default_na=False).drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S16_Zero_fit", index=False)
    pd.read_csv(R / "zero_shape_factorial_summary.csv",keep_default_na=False).drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S16_Zero_prior_factorial", index=False)
    pd.read_csv(R / "publication_DV_comparators.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S17_DV_comparators", index=False)
    covariance_diagnosis.assign(analysis_role="gene-wise covariance diagnostic").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
        writer, sheet_name="S18_Covariance", index=False)
    pd.read_csv(R / "current_simulation_by_abundance.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S3_Abundance", index=False)
    pd.read_csv(R / "CV2_estimator_comparison.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S4_CV2_estimators", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_distribution_stability_summary.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S15_Distribution_stability", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_stability_coverage.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S15_Mapping", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_example_features.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S15_Examples", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_enrichment_summary.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S21_Enrichment_summary", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_enrichment_mapping.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S21_Annotation_mapping", index=False)
    pd.read_csv(R / f"sgcb_{VERSION}_kirc_enrichment.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S21_All_tested_terms", index=False)
    pd.read_csv(R / "publication_zero_calibration.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S22_Zero_calibration", index=False)
    pd.read_csv(R / "publication_zero_calibration_scenarios.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S22_Zero_scenarios", index=False)
    pd.read_csv(R / "zero_calibration_parameters.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S22_Parameters", index=False)
    pd.read_csv(R / "zero_calibration_truth.csv").drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S22_Truth", index=False)
    new_biology.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S24_Disease_cohorts",index=False)
    new_covariates.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S24_Covariate_models",index=False)
    pd.concat([
        pd.read_csv(ROOT/"data/new_biology/sle_classical_monocytes_2021/samples.csv").assign(dataset=COHORTS[0]),
        pd.read_csv(R/"current_masld_liver_2023_covariate_samples.csv").assign(dataset=COHORTS[1])],
        ignore_index=True).drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S24_Samples",index=False)
    pd.read_csv(R/"current_masld_liver_2023_composition_markers.csv").assign(
        dataset=COHORTS[1]).drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S24_Markers",index=False)
    pd.concat(map(pd.read_csv,bio_prefix+"_example_summary.csv"),keys=COHORTS,
        names=["source_dataset","record"]).reset_index().drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
            writer,sheet_name="S24_Examples",index=False)
    pd.concat(map(pd.read_csv,bio_prefix+"_example_effects.csv"),keys=COHORTS,
        names=["source_dataset","record"]).reset_index().drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
            writer,sheet_name="S24_Example_effects",index=False)
    pd.concat(map(pd.read_csv,bio_prefix+"_example_empirical.csv"),keys=COHORTS,
        names=["source_dataset","record"]).reset_index().drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
            writer,sheet_name="S24_Example_moments",index=False)
    pd.concat(map(pd.read_csv,bio_prefix+"_enrichment_summary.csv"),keys=COHORTS,
        names=["source_dataset","record"]).reset_index().drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
            writer,sheet_name="S25_Enrichment_summary",index=False)
    pd.concat(map(pd.read_csv,bio_prefix+"_enrichment_mapping.csv"),keys=COHORTS,
        names=["source_dataset","record"]).reset_index().drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(
            writer,sheet_name="S25_Annotation_mapping",index=False)
    new_enrichment.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer,sheet_name="S25_All_tested_terms",index=False)
    population_calibration.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Population_calibration", index=False)
    population_scenarios.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Population_scenarios", index=False)
    population_family_summary.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Kernel_families", index=False)
    population_family_scenarios.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Kernel_scenarios", index=False)
    population_family_operating.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Kernel_rank_curves", index=False)
    population_biology.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Population_biology", index=False)
    population_covariates.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Population_covariates", index=False)
    population_biology_genes.drop(columns=CONFIG["publication_metadata_columns"], errors="ignore").replace(CONFIG["publication_label_replacements"], regex=True).to_excel(writer, sheet_name="S26_Population_genes", index=False)
print(null_summary.to_string(index=False))
print(table2.to_string(index=False))
print(table4.to_string(index=False))
