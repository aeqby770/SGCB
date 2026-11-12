from pathlib import Path
import os
import sys
import runpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.graphics.mosaicplot import mosaic
from matplotlib.patches import ConnectionPatch, FancyBboxPatch, FancyArrowPatch, Rectangle, Patch
from matplotlib.colors import LogNorm, SymLogNorm, LinearSegmentedColormap, TwoSlopeNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
from textwrap import fill
from scipy.special import gammaln
from scipy.stats import t, spearmanr
from sklearn.metrics import precision_recall_curve

CONFIG = {"paths": {"root": Path(sys.argv[1]),
    "previous_results": Path(sys.argv[1]) / "tmp/revision_0.5.5_integrated/results"}, "threads": int(os.environ["SGCB_THREADS"]), "plot": {
    "font.family": ["Arial", "DejaVu Sans"], "font.size": 9.3, "axes.titlesize": 9.5,
    "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.titlepad": 9,
    "axes.labelsize": 9.2, "xtick.labelsize": 8.2, "ytick.labelsize": 8.2, "legend.fontsize": 8,
    "text.color": "#192E38", "axes.labelcolor": "#192E38",
    "xtick.color": "#192E38", "ytick.color": "#192E38",
    "axes.edgecolor": "#627984", "axes.axisbelow": True,
    "grid.color": "#E8EDEE", "grid.linewidth": .45,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": .55, "ytick.major.width": .55,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "mathtext.fontset": "custom", "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic", "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial", "mathtext.fallback": "stix",
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
    "lines.linewidth": 1.55, "lines.markersize": 5.5, "lines.markeredgewidth": .6,
    "lines.solid_capstyle": "round", "errorbar.capsize": 2.7,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 300, "figure.constrained_layout.h_pad": 4 / 72,
    "figure.constrained_layout.w_pad": 4 / 72,
    "figure.constrained_layout.hspace": 0, "figure.constrained_layout.wspace": 0},
    "figure_palettes": {'1': {'name': 'Lake, coral and cream',
           'colors': {'DE': '#3B9DB5',
                      'DV': '#E98673',
                      'both': '#E8D6A5',
                      'normal': '#3B9DB5',
                      'tumor': '#E98673',
                      'reference': '#BDCECE',
                      'ink': '#253B45',
                      'false': '#AD5260',
                      'other': '#D9E5E1'},
           'cohorts': {'normal': '#ACDCE4', 'tumor': '#F3B5A4', 'heldout_tumor': '#F7EACB'},
           'light': '#FFF8E9',
           'soft': '#BADFE3'},
     '2': {'name': 'Slate blue and amber',
           'colors': {'DE': '#6E95AE',
                      'DV': '#3D7194',
                      'both': '#CFAB68',
                      'normal': '#8EA3B1',
                      'tumor': '#3D7194',
                      'reference': '#CDD7DD',
                      'ink': '#253B45',
                      'false': '#AD5260',
                      'other': '#E7EDF0'},
           'cohorts': {'normal': '#B8CCD7', 'tumor': '#7CA5BC', 'heldout_tumor': '#E8D6AA'},
           'light': '#FBF7EC',
           'soft': '#B4CFDD',
           'bias_scale': ['#437A9A', '#C8DDE8', '#FCF9EF', '#EDDBB1', '#C69C55'],
           'error_scale': ['#FFF7DF', '#EBCB89', '#DA9676', '#BF6D57', '#923F3D']},
     '3': {'name': 'Sage and apricot',
           'colors': {'DE': '#7CA794',
                      'DV': '#C7926B',
                      'both': '#D7BF8A',
                      'normal': '#7CA794',
                      'tumor': '#C7926B',
                      'reference': '#D5DDD6',
                      'ink': '#253B45',
                      'false': '#B7747C',
                      'other': '#E8ECE5',
                      'control': '#7F948C',
                      'case': '#C7926B'},
           'cohorts': {'normal': '#C8DFD0', 'tumor': '#E9C3A6', 'heldout_tumor': '#EFE5CB'},
           'light': '#FBFAF1',
           'soft': '#C6DCCB',
           'availability': ['#FFF5D8', '#D9DEB7', '#8FC6B2', '#4B9D9D', '#286677'],
           'calls': ['#EDF5F2', '#F3E0AE', '#E8BA80', '#D78A6D', '#A44E45']},
     '4': {'name': 'Pomegranate and mist blue',
           'colors': {'DE': '#AD5972',
                      'DV': '#7198B4',
                      'both': '#CFB27A',
                      'normal': '#7198B4',
                      'tumor': '#AD5972',
                      'reference': '#C8D5DD',
                      'ink': '#253B45',
                      'false': '#AD5260',
                      'other': '#E8EDF0'},
           'cohorts': {'normal': '#C9DCE8', 'tumor': '#E9C4CE', 'heldout_tumor': '#EFE1BE'},
           'light': '#FCF5F7',
           'soft': '#E6C5CF',
           'methods': {'SGCB': '#AD5972',
                       'DESeq2': '#628EAD',
                       'edgeR': '#83A88E',
                       'limma': '#C8AA6C',
                       'glmGamPoi': '#9C8C79',
                       'NOISeq': '#A8BEC5',
                       'EBSeq': '#D1A6AE'}},
     '5': {'name': 'Turquoise and sand',
           'colors': {'DE': '#809BA9',
                      'DV': '#4C9F98',
                      'both': '#D5B77B',
                      'normal': '#8EBDB1',
                      'tumor': '#4C9F98',
                      'reference': '#C9DAD5',
                      'ink': '#253B45',
                      'false': '#AD5260',
                      'other': '#E7EFEA'},
           'cohorts': {'normal': '#C4DFD5', 'tumor': '#8AC2B8', 'heldout_tumor': '#EBDCB6'},
           'light': '#F9F8EF',
           'soft': '#C3DED3',
           'density_scale': ['#EAF3EF', '#B4D5C8', '#72B1A0', '#2E8078']},
     '6': {'name': 'Sky blue and rose',
           'colors': {'DE': '#72A8C1',
                      'DV': '#CB829B',
                      'both': '#D9C39A',
                      'normal': '#93B5CC',
                      'tumor': '#CB829B',
                      'reference': '#CBD8E1',
                      'ink': '#253B45',
                      'false': '#AD5260',
                      'other': '#DCE7EB',
                      'control': '#93B5CC',
                      'case': '#CB829B'},
           'cohorts': {'normal': '#C5D9E8', 'tumor': '#E9BFCE', 'heldout_tumor': '#EFE1C7'},
           'light': '#FCF5F8',
           'soft': '#D6E4ED'}},
    "hexbin": {"gridsize": 36, "mincnt": 1, "linewidths": 0, "rasterized": False},
    "designs": ["10:10", "10:50", "50:50"], "alpha": .05,
    "mean_record_id": os.environ.get("SGCB_MEAN_RECORD_ID", "0.5.7"),
    "distribution": {"bandwidth_adjust": .85, "width": .72, "point_size": .70},
    "kirc_enrichment_BP_terms": 3,
    "kirc_enrichment_terms": ["GO:0008380", "GO:0032543", "GO:0043161", "GO:0005761"], "kirc_enrichment_label_width": 30,
    "kirc_enrichment_labels": {"GO:0008380": "RNA splicing", "GO:0032543": "Mitochondrial\ntranslation", "GO:0043161": "Proteasome-mediated\nprotein degradation", "GO:0005761": "Mitochondrial ribosome"},
    "family_comparison": {"figsize": [7.25, 6.9], "height_ratios": [.55, 1.10, 1.40],
                          "curve_widths": [1.5, .75, .75]},
    "heatmap": {"cell_gap": 2.1, "font_size": 8, "colorbar_height": .035,
                "header_bounds": [0, 1.012, 1, .022], "bias_linear_range": .05},
    "kirc_zoom": {"bounds": [.05,.06,.42,.32], "quantiles": [.025,.975]},
    "detail": {"ink": "#192E38", "outline": "#4D6975", "point_edge": "white",
               "point_edge_width": .45, "interval_width": 1.5, "interval_caps": 2.7,
               "bar_edge_width": .45, "legend_size": 7.8},
    "resource_datasets": ["channel_n10_10","channel_n10_50","channel_n50_50","zero_n10_10","zero_n10_50","zero_n50_50","readcount_null_1","readcount_null_2","readcount_null_3_streamed","readcount_null_4","readcount_null_5","readcount_null_6","readcount_null_7","readcount_null_8","readcount_null_9","readcount_null_10","sle_classical_monocytes_2021","kirc","masld_liver_2023"],
    "resource_cohorts": ["sle_classical_monocytes_2021","kirc","masld_liver_2023"],
    "resource_marker": {"offset":8,"per_sample":.3},
    "resource_cohort_labels": ["SLE\nmonocytes","KIRC","MASLD\nliver"]}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
plt.rcParams.update(CONFIG["plot"])
ROOT = CONFIG["paths"]["root"]
RECORD_ID = os.environ.get("SGCB_STUDY_RECORD_ID", "0.5.9")
R, F = ROOT / "results", ROOT / "manuscript" / "figures"
F.mkdir(parents=True, exist_ok=True)
INK=CONFIG["detail"]["ink"]
RED=CONFIG["figure_palettes"]["1"]["colors"]["false"]
CHANNELS = ["DE_only", "DV_only", "DE_and_DV", "neither"]
CHANNEL_LABELS = {"DE_only": "DE only", "DV_only": "DV only", "DE_and_DV": "Both", "neither": "Neither"}
SCENARIOS = ["null_gg", "mean", "variance", "mean_variance", "shape", "null_lognormal",
             "null_contamination", "null_mixture", "shape_mixture"]
SCENARIO_LABELS = ["GG null", "Mean shift", "CV² shift", "Mean + CV²", "Equal-moment shape",
                  "Lognormal null", "Contaminated null", "Mixture null", "Equal-moment mixture"]


FIG_STYLE=CONFIG["figure_palettes"]["1"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure 1 presents the mechanism as connected input, model and readout objects.
# The Figure 1 renderer uses a local font and canvas style.
with plt.rc_context():
    runpy.run_path(str(ROOT / "visualization" / "overview.py"),
        init_globals={"ROOT":ROOT,"FIG_STYLE":FIG_STYLE,"THREADS":CONFIG["threads"]})

# Kernel restrictions share the GG fits and fixed input designs.
population_family_inputs = pd.MultiIndex.from_product(
    [["GG","Gamma","Lognormal"],["channel_simulation","zero_calibration"],
     ["10_10","10_50","50_50"]],
    names=["family","input_panel","sample_design"]).to_frame(index=False)
population_family_inputs["fit"] = population_family_inputs.index
population_family_inputs["panel"] = population_family_inputs.input_panel.map(
    {"channel_simulation":"Positive","zero_calibration":"Observed zeros"})
population_family_inputs["prefix"] = np.where(
    population_family_inputs.family.eq("GG"),
    population_family_inputs.input_panel + "_population_mixture_n" +
        population_family_inputs.sample_design + "_refined_continued_paired",
    "population_" + population_family_inputs.family.str.lower() + "_" +
        population_family_inputs.input_panel + "_n" + population_family_inputs.sample_design)
population_family_inputs["prefix"] = str(R) + "/" + population_family_inputs.prefix
population_family_metadata = population_family_inputs[["fit","family","panel"]]
# The public API supplies the arithmetic-mean and CV-squared endpoints.
# The archived fits below provide historical schema and separate sensitivity records.
current_inputs = population_family_inputs.loc[
    population_family_inputs.family.eq("GG"), ["panel", "sample_design"]].copy()
current_inputs["prefix"] = (str(R) + f"/sgcb_{RECORD_ID}_" +
    current_inputs.panel.map({"Positive": "channel", "Observed zeros": "zero"}) +
    "_n" + current_inputs.sample_design)
api_metrics = pd.concat(map(pd.read_csv, current_inputs.prefix + "_metrics.csv"),
    ignore_index=True)
api_scenarios = pd.concat(map(pd.read_csv, current_inputs.prefix + "_scenario_metrics.csv"),
    ignore_index=True)
api_genes = pd.concat(map(pd.read_csv, current_inputs.prefix + "_evaluation_genes.csv"),
    ignore_index=True)
assert len(api_metrics) == 12 and api_metrics[["panel", "design"]].drop_duplicates().shape[0] == 6
api_metrics.to_csv(R / "publication_current_api_metrics.csv", index=False, encoding="utf-8")
api_scenarios.to_csv(R / "publication_current_api_scenarios.csv", index=False, encoding="utf-8")

final_summary = api_metrics.loc[api_metrics.endpoint.eq("DV")].copy()
final_summary["testing_family_size"] = final_summary.features
final_summary["features"] = final_summary.finite_truth
final_summary["available"] = final_summary.interval_available
final_summary["availability"] = final_summary.available / final_summary.features
final_summary["calls"] = final_summary.TP + final_summary.FP + final_summary.undefined_calls
final_summary["raw_null_rejections"] = np.rint(final_summary.raw_null_rate * final_summary.nulls).astype(int)
final_summary["coverage_interval_enclosure"] = final_summary.coverage_available
final_summary["coverage_available"] = final_summary.native_coverage_available
final_summary["coverage_all"] = final_summary.native_coverage_available * final_summary.availability
final_summary = final_summary.rename(columns={"nulls": "null_features",
    "alternatives": "alternative_features", "undefined_calls": "calls_with_undefined_truth",
    "raw_null_rate": "raw_null_rate_all"})
final_summary["truth_domain"] = "finite_CV2"
final_summary["method"] = "EB-GG-paired-population"
final_summary["family"] = "GG"
final_summary["panel"] = final_summary.panel.map({"channel": "Positive", "zero": "Observed zeros"})
final_summary["sensitivity"] = final_summary.power
final_scenarios = api_scenarios.loc[api_scenarios.endpoint.eq("DV")].copy()
final_scenarios["truth_domain"] = np.where(final_scenarios.truth_finite.gt(0),
    "finite_CV2", "undefined_CV2")
final_scenarios["coverage_interval_enclosure"] = final_scenarios.coverage_available
final_scenarios["coverage_available"] = final_scenarios.native_coverage_available
final_scenarios["coverage_all"] = final_scenarios.native_coverage_available * (
    final_scenarios.available / final_scenarios.features)
final_scenarios["raw_p005_available"] = final_scenarios.raw_p005 * (
    final_scenarios.features / final_scenarios.available)
final_scenarios = final_scenarios.rename(columns={"raw_p005": "raw_p005_all"})
final_scenarios["method"] = "EB-GG-paired-population"
final_scenarios["family"] = "GG"
final_scenarios["panel"] = final_scenarios.panel.map({"channel": "Positive", "zero": "Observed zeros"})
final_genes = api_genes.loc[api_genes.endpoint.eq("DV")].copy().rename(
    columns={"pvalue": "p", "qvalue": "q", "truth": "dv_truth",
             "covered": "interval_enclosure_covered", "native_covered": "covered"})
final_genes["method"] = "EB-GG-paired-population"
final_genes["family"] = "GG"
final_genes["panel"] = final_genes.panel.map({"channel": "Positive", "zero": "Observed zeros"})
# Kernel-family comparisons use matched fits from the same package version.
matched_prefixes = (str(R) + f"/cv_matched_families_{RECORD_ID}_" +
    current_inputs.panel.map({"Positive": "channel", "Observed zeros": "zero"}) +
    "_n" + current_inputs.sample_design)
matched_metadata = current_inputs[["panel", "sample_design"]].copy()
matched_metadata["fit"] = current_inputs.panel + "|" + current_inputs.sample_design
matched_metadata["design"] = current_inputs.sample_design.str.replace("_", ":", regex=False)
matched_genes = pd.concat(map(pd.read_csv, matched_prefixes + "_genes.csv"),
    keys=matched_metadata.fit, names=["fit", "record"]).reset_index().merge(
        matched_metadata[["fit", "panel", "design"]], on="fit", validate="many_to_one")
matched_genes["family"] = matched_genes.method.replace({"lognormal": "Lognormal"})
matched_genes = matched_genes.rename(columns={"pvalue": "p", "adjusted": "q", "truth": "dv_truth"})
matched_genes["valid"] = matched_genes.available & matched_genes.truth_finite
matched_genes["is_null"] = matched_genes.truth_finite & ~matched_genes.dv_truth
matched_genes["is_alternative"] = matched_genes.truth_finite & matched_genes.dv_truth
matched_genes["true_call"] = matched_genes.called & matched_genes.is_alternative
matched_genes["false_call"] = matched_genes.called & matched_genes.is_null
matched_genes["undefined_call"] = matched_genes.called & ~matched_genes.truth_finite
matched_genes["raw_null_call"] = matched_genes.is_null & matched_genes.p.lt(CONFIG["alpha"])
matched_genes["raw_call"] = matched_genes.p.lt(CONFIG["alpha"])
matched_genes["valid_null"] = matched_genes.valid & matched_genes.is_null
matched_genes["valid_covered"] = matched_genes.valid & matched_genes.covered
matched_genes["error_valid"] = np.where(matched_genes.valid,
    matched_genes.effect - matched_genes.truth_effect, np.nan)
matched_genes["squared_error_valid"] = matched_genes.error_valid ** 2
matched_genes["se_valid"] = np.where(matched_genes.valid, matched_genes.se, np.nan)
family_aggregations = {
    "testing_family_size": ("gene_id", "size"), "finite_truth": ("truth_finite", "sum"),
    "available_all": ("available", "sum"), "available_finite": ("valid", "sum"),
    "null_features": ("is_null", "sum"), "alternative_features": ("is_alternative", "sum"),
    "TP": ("true_call", "sum"), "FP": ("false_call", "sum"),
    "calls_with_undefined_truth": ("undefined_call", "sum"),
    "raw_null_rejections": ("raw_null_call", "sum"), "raw_rejections": ("raw_call", "sum"),
    "valid_nulls": ("valid_null", "sum"), "covered_count": ("valid_covered", "sum"),
    "RMSE": ("squared_error_valid", "mean"), "bias": ("error_valid", "mean"),
    "median_se": ("se_valid", "median")}
family_overall = matched_genes.groupby(["panel", "design", "family"], sort=False).agg(
    **family_aggregations).reset_index().assign(scope="overall")
family_by_scenario = matched_genes.groupby(["panel", "design", "family", "scenario"],
    sort=False).agg(**family_aggregations).reset_index().assign(scope="scenario")
family_results = pd.concat([family_overall, family_by_scenario], ignore_index=True)
family_results["features"] = np.where(family_results.scope.eq("overall"),
    family_results.finite_truth, family_results.testing_family_size)
family_results["available"] = np.where(family_results.scope.eq("overall"),
    family_results.available_finite, family_results.available_all)
family_results["RMSE"] = np.sqrt(family_results.RMSE)
family_results["calls"] = family_results.TP + family_results.FP + family_results.calls_with_undefined_truth
family_results["power"] = family_results.TP / family_results.alternative_features
family_results["sensitivity"] = family_results.power
family_results["FDP"] = family_results.FP / np.maximum(1, family_results.TP + family_results.FP)
family_results["availability"] = family_results.available / family_results.features
family_results["coverage_available"] = family_results.covered_count / family_results.available_finite
family_results["coverage_all"] = family_results.covered_count / family_results.finite_truth
family_results["raw_null_rate_all"] = family_results.raw_null_rejections / family_results.null_features
family_results["raw_null_rate_available"] = family_results.raw_null_rejections / family_results.valid_nulls
family_results["raw_p005_all"] = family_results.raw_rejections / family_results.testing_family_size
family_results["raw_p005_available"] = family_results.raw_rejections / family_results.available_all
family_results["truth_domain"] = np.where(family_results.finite_truth.gt(0), "finite_CV2", "undefined_CV2")
family_results["method"] = "EB-" + family_results.family + "-paired-population"
matched_times = pd.concat(map(pd.read_csv,matched_prefixes + "_metrics.csv"),
    keys=matched_metadata.fit,names=["fit","record"]).reset_index().merge(
        matched_metadata[["fit","panel","design"]],on="fit",validate="many_to_one")
matched_times["family"] = matched_times.method.replace({"lognormal":"Lognormal"})
population_family_summary = family_results.loc[family_results.scope.eq("overall")].copy().merge(
    matched_times[["panel","design","family","seconds","timing_source"]],
    on=["panel","design","family"],validate="one_to_one")
population_family_scenarios = family_results.loc[family_results.scope.eq("scenario")].copy()
population_family_genes = matched_genes.copy()


# Figure 2 ranks all available positive-panel scores; the power denominator stays 1200.
population_ranked = population_family_genes.loc[
    population_family_genes.panel.eq("Positive")].copy()
population_ranked["rank_TP"] = population_ranked.dv_truth.astype(int)
population_ranked["rank_FP"] = (~population_ranked.dv_truth.astype(bool)).astype(int)
population_family_operating = population_ranked.groupby(["design","family","p"],sort=True).agg(
    TP=("rank_TP","sum"),FP=("rank_FP","sum")).reset_index()
population_family_operating[["TP","FP"]] = population_family_operating.groupby(
    ["design","family"],sort=False)[["TP","FP"]].cumsum()
population_family_operating = population_family_operating.merge(
    population_family_summary.loc[population_family_summary.panel.eq("Positive"),
        ["design","family","alternative_features"]],on=["design","family"],validate="many_to_one")
population_family_operating["sensitivity"] = (
    population_family_operating.TP / population_family_operating.alternative_features)
population_family_operating["FDP"] = population_family_operating.FP / (
    population_family_operating.TP + population_family_operating.FP)
# Discrimination uses all finite-truth features; unavailable scores rank last.
api_de_rank = api_genes.loc[api_genes.endpoint.eq("DE"),
    ["panel","design","endpoint","gene_id","truth","truth_finite","pvalue"]].rename(columns={"pvalue":"p"}).copy()
api_de_rank["panel"] = api_de_rank.panel.map({"channel":"Positive","zero":"Observed zeros"})
api_de_rank["family"] = "GG"
family_dv_rank = matched_genes[["panel","design","family","gene_id","dv_truth","truth_finite","p"]].rename(
    columns={"dv_truth":"truth"}).assign(endpoint="DV")
rank_data = pd.concat([api_de_rank,family_dv_rank],ignore_index=True)
rank_data = rank_data.loc[rank_data.truth_finite].copy()
rank_data["ranking_p"] = rank_data.p.fillna(1)
rank_keys = ["panel","design","endpoint","family"]
rank_bins = rank_data.groupby(rank_keys+["ranking_p"],sort=True).agg(
    true_at_tie=("truth","sum"),features_at_tie=("gene_id","size")).reset_index()
rank_bins["cum_true"] = rank_bins.groupby(rank_keys,sort=False).true_at_tie.cumsum()
rank_bins["cum_features"] = rank_bins.groupby(rank_keys,sort=False).features_at_tie.cumsum()
rank_bins["true_total"] = rank_bins.groupby(rank_keys,sort=False).true_at_tie.transform("sum")
rank_bins["AP_increment"] = (rank_bins.cum_true/rank_bins.cum_features *
    rank_bins.true_at_tie/rank_bins.true_total)
discrimination = rank_bins.groupby(rank_keys,sort=False).AP_increment.sum().rename(
    "average_precision_finite_truth").reset_index()
population_family_summary = population_family_summary.merge(
    discrimination.loc[discrimination.endpoint.eq("DV")].drop(columns="endpoint"),
    on=["panel","design","family"],validate="one_to_one")
population_family_summary["precision_finite_truth"] = (
    population_family_summary.TP/(population_family_summary.TP+population_family_summary.FP).clip(lower=1))
population_family_summary["specificity_finite_truth"] = (
    1-population_family_summary.FP/population_family_summary.null_features)
population_family_summary["F1_finite_truth"] = 2*population_family_summary.TP/(
    population_family_summary.alternative_features+population_family_summary.TP+
    population_family_summary.FP).clip(lower=1)
api_discrimination = discrimination.loc[discrimination.family.eq("GG")].drop(columns="family").copy()
api_discrimination["panel"] = api_discrimination.panel.map({"Positive":"channel","Observed zeros":"zero"})
api_metrics = api_metrics.merge(api_discrimination,on=["panel","design","endpoint"],validate="one_to_one")
api_metrics["precision_finite_truth"] = api_metrics.TP/(api_metrics.TP+api_metrics.FP).clip(lower=1)
api_metrics["specificity_finite_truth"] = 1-api_metrics.FP/api_metrics.nulls
api_metrics["F1_finite_truth"] = 2*api_metrics.TP/(
    api_metrics.alternatives+api_metrics.TP+api_metrics.FP).clip(lower=1)
api_metrics.to_csv(R / "publication_current_api_metrics.csv",index=False,encoding="utf-8")

population_family_summary.to_csv(R / "publication_population_family_summary.csv",
    index=False,encoding="utf-8")
population_family_scenarios.to_csv(R / "publication_population_family_scenarios.csv",
    index=False,encoding="utf-8")
population_family_operating.to_csv(R / "publication_population_family_operating_curves.csv",
    index=False,encoding="utf-8")


FIG_STYLE=CONFIG["figure_palettes"]["2"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure 2: matched columns separate the fixed decision from truth-dependent ranking.
family_summary = population_family_summary.loc[population_family_summary.panel.eq("Positive")]
family_scenarios = population_family_scenarios.loc[population_family_scenarios.panel.eq("Positive")]
operating = population_family_operating
families = ["Gamma", "Lognormal", "GG"]
family_colors = [NORMAL, JOINT_COLOR, DV_COLOR]
design_order = CONFIG["designs"]
s0 = family_summary.query("design == '10:10'").set_index("family").loc[families]
s1 = family_summary.query("design == '10:50'").set_index("family").loc[families]
s2 = family_summary.query("design == '50:50'").set_index("family").loc[families]

fig = plt.figure(figsize=CONFIG["family_comparison"]["figsize"], layout="constrained")
top, middle, bottom = fig.subfigures(3, 1, height_ratios=CONFIG["family_comparison"]["height_ratios"])
top.suptitle("$\\mathbf{a}$",x=.015,ha="left",fontsize=10)
a0,a1,a2=top.subplots(1,3,sharex=True,sharey=True)
a0.scatter(s0.TP,[2,1,0],c=family_colors,s=28,zorder=3)
a1.scatter(s1.TP,[2,1,0],c=family_colors,s=28,zorder=3)
a2.scatter(s2.TP,[2,1,0],c=family_colors,s=28,zorder=3)
np.frompyfunc(a0.text,3,1)(np.repeat(1350,3),[2,1,0],np.char.mod("%d",s0.FP))
np.frompyfunc(a1.text,3,1)(np.repeat(1350,3),[2,1,0],np.char.mod("%d",s1.FP))
np.frompyfunc(a2.text,3,1)(np.repeat(1350,3),[2,1,0],np.char.mod("%d",s2.FP))
plt.setp(list(a0.texts)+list(a1.texts)+list(a2.texts),ha="center",va="center",fontsize=8,color=RED)
plt.setp([a0,a1,a2],xlim=(0,1450),ylim=(-.6,2.6),yticks=[2,1,0],xticks=[0,600,1200],xlabel="True calls")
a0.set(title="10 : 10",yticklabels=["Gamma","LN","GG"])
a1.set(title="10 : 50")
a2.set(title="50 : 50")
a0.text(1350,2.8,"FP",ha="center",fontsize=7.5,color=RED)
a1.text(1350,2.8,"FP",ha="center",fontsize=7.5,color=RED)
a2.text(1350,2.8,"FP",ha="center",fontsize=7.5,color=RED)
plt.setp([a0,a1,a2],axisbelow=True)
a0.grid(axis="x",color="#E6ECEF",lw=.6)
a1.grid(axis="x",color="#E6ECEF",lw=.6)
a2.grid(axis="x",color="#E6ECEF",lw=.6)

middle.suptitle("$\\mathbf{b}$", x=.015, ha="left", fontsize=10, weight="normal")
b0, b1, b2 = middle.subplots(1, 3, sharey=True,
    gridspec_kw={"width_ratios":CONFIG["family_comparison"]["curve_widths"]})
d0 = operating.query("design == '10:10'")
d1 = operating.query("design == '10:50'")
d2 = operating.query("design == '50:50'")
b0.plot(d0.loc[d0.family.eq("Gamma"), "FDP"], d0.loc[d0.family.eq("Gamma"), "sensitivity"], color=NORMAL)
b0.plot(d0.loc[d0.family.eq("Lognormal"), "FDP"], d0.loc[d0.family.eq("Lognormal"), "sensitivity"], color=JOINT_COLOR)
b0.plot(d0.loc[d0.family.eq("GG"), "FDP"], d0.loc[d0.family.eq("GG"), "sensitivity"], color=DV_COLOR, lw=1.9)
b1.plot(d1.loc[d1.family.eq("Gamma"), "FDP"], d1.loc[d1.family.eq("Gamma"), "sensitivity"], color=NORMAL)
b1.plot(d1.loc[d1.family.eq("Lognormal"), "FDP"], d1.loc[d1.family.eq("Lognormal"), "sensitivity"], color=JOINT_COLOR)
b1.plot(d1.loc[d1.family.eq("GG"), "FDP"], d1.loc[d1.family.eq("GG"), "sensitivity"], color=DV_COLOR, lw=1.9)
b2.plot(d2.loc[d2.family.eq("Gamma"), "FDP"], d2.loc[d2.family.eq("Gamma"), "sensitivity"], color=NORMAL)
b2.plot(d2.loc[d2.family.eq("Lognormal"), "FDP"], d2.loc[d2.family.eq("Lognormal"), "sensitivity"], color=JOINT_COLOR)
b2.plot(d2.loc[d2.family.eq("GG"), "FDP"], d2.loc[d2.family.eq("GG"), "sensitivity"], color=DV_COLOR, lw=1.9)
b0.plot(s0.loc["Gamma", "FDP"], s0.loc["Gamma", "sensitivity"], "s", color=NORMAL, ms=5, zorder=6)
b0.plot(s0.loc["Lognormal", "FDP"], s0.loc["Lognormal", "sensitivity"], "^", color=JOINT_COLOR, mfc="white", ms=8, zorder=6)
b0.plot(s0.loc["GG", "FDP"], s0.loc["GG", "sensitivity"], "o", color=DV_COLOR, ms=3.5, zorder=7)
b1.plot(s1.loc["Gamma", "FDP"], s1.loc["Gamma", "sensitivity"], "s", color=NORMAL, ms=5, zorder=6)
b1.plot(s1.loc["Lognormal", "FDP"], s1.loc["Lognormal", "sensitivity"], "^", color=JOINT_COLOR, mfc="white", ms=6, zorder=6)
b1.plot(s1.loc["GG", "FDP"], s1.loc["GG", "sensitivity"], "o", color=DV_COLOR, ms=4.5, zorder=7)
b2.plot(s2.loc["Gamma", "FDP"], s2.loc["Gamma", "sensitivity"], "s", color=NORMAL, ms=5, zorder=6)
b2.plot(s2.loc["Lognormal", "FDP"], s2.loc["Lognormal", "sensitivity"], "^", color=JOINT_COLOR, mfc="white", ms=6, zorder=6)
b2.plot(s2.loc["GG", "FDP"], s2.loc["GG", "sensitivity"], "o", color=DV_COLOR, ms=4.5, zorder=7)
b0.grid(True)
b1.grid(True)
b2.grid(True)
b0.axvline(.05, color=RED, ls=":", lw=.8)
b1.axvline(.05, color=RED, ls=":", lw=.8)
b2.axvline(.05, color=RED, ls=":", lw=.8)
b0.set_xscale("symlog", linthresh=.05, linscale=.5)
b1.set_xscale("symlog", linthresh=.05, linscale=.5)
b2.set_xscale("symlog", linthresh=.05, linscale=.5)
plt.setp([b0, b1, b2], xlim=(-.005, 1.03), ylim=(-.04, 1.04),
         xticks=[0, .05, .2, 1], xticklabels=["0", ".05", ".2", "1"], xlabel="Realized FDP")
b0.set(title="10 : 10", ylabel="Sensitivity")
b1.set(title="10 : 50")
b2.set(title="50 : 50")
middle.legend(handles=[Line2D([], [], color=NORMAL, marker="s", ms=4, label="Gamma"),
    Line2D([], [], color=JOINT_COLOR, marker="^", mfc="white", ms=5, label="Lognormal"),
    Line2D([], [], color=DV_COLOR, marker="o", ms=4, label="GG")],
    title="Markers: default BH", title_fontsize=7.3, frameon=False,
    loc="outside lower center", ncol=3, fontsize=7.6, handlelength=1.7)

error_grid=bottom.add_gridspec(2,2,height_ratios=[1,CONFIG["heatmap"]["colorbar_height"]])
c=bottom.add_subplot(error_grid[0,0])
d=bottom.add_subplot(error_grid[0,1],sharey=c)
bias_bar=bottom.add_subplot(error_grid[1,0])
rmse_bar=bottom.add_subplot(error_grid[1,1])
bias50=family_scenarios.query("design == '50:50'").pivot(index="scenario",columns="family",values="bias").reindex(SCENARIOS)[families]
rmse50=family_scenarios.query("design == '50:50'").pivot(index="scenario",columns="family",values="RMSE").reindex(SCENARIOS)[families]
bias_limit=np.abs(bias50.to_numpy()).max()
bias_ticks=np.array([-bias_limit,-.1,0,.1,bias_limit])
rmse_limits=[float(rmse50.to_numpy().min()),float(rmse50.to_numpy().max())]
rmse_ticks=np.geomspace(*rmse_limits,3)
bias_labels=np.where(np.abs(bias50.to_numpy())<.005,"0",np.char.mod("%+.2f",bias50.to_numpy()))
bias_palette=sns.blend_palette(FIG_STYLE["bias_scale"],as_cmap=True)
error_palette=sns.blend_palette(FIG_STYLE["error_scale"],as_cmap=True)
sns.heatmap(bias50,ax=c,annot=bias_labels,fmt="",cmap=bias_palette,
    norm=SymLogNorm(linthresh=CONFIG["heatmap"]["bias_linear_range"],vmin=-bias_limit,vmax=bias_limit),
    linewidths=CONFIG["heatmap"]["cell_gap"],linecolor="white",
    cbar_ax=bias_bar,cbar_kws={"orientation":"horizontal","ticks":bias_ticks,"label":"Log₂-effect bias (50:50)"},
    annot_kws={"fontsize":CONFIG["heatmap"]["font_size"]},xticklabels=["Gamma","LN","GG"],yticklabels=False)
sns.heatmap(rmse50,ax=d,annot=True,fmt=".2f",cmap=error_palette,
    norm=LogNorm(vmin=rmse_limits[0],vmax=rmse_limits[1]),linewidths=CONFIG["heatmap"]["cell_gap"],linecolor="white",
    cbar_ax=rmse_bar,cbar_kws={"orientation":"horizontal","ticks":rmse_ticks,"label":"Log₂-effect RMSE (50:50)"},
    annot_kws={"fontsize":CONFIG["heatmap"]["font_size"]},xticklabels=["Gamma","LN","GG"],yticklabels=False)
c.set(title=r"$\mathbf{c}$",xlabel="",ylabel="",yticks=np.arange(9)+.5,
    yticklabels=["Gamma null","Mean shift (Gamma)","CV² shift (Gamma)","Mean + CV² (Gamma)",
    "GG → GG: equal CV²","Lognormal null","Contaminated null","Mixture null","GG → mixture: equal CV²"])
d.set(title=r"$\mathbf{d}$",xlabel="",ylabel="")
c.tick_params(length=0,labelrotation=0)
d.tick_params(length=0,labelrotation=0,labelleft=False)
bias_bar.set_xticklabels(np.char.mod("%.2g",bias_ticks))
rmse_bar.set_xticklabels(np.char.mod("%.2g",rmse_ticks))
bias_bar.tick_params(labelsize=7,length=2,pad=3)
rmse_bar.tick_params(labelsize=7,length=2,pad=3)
bias_bar.minorticks_off()
rmse_bar.minorticks_off()
bias_bar.xaxis.label.set_size(8)
rmse_bar.xaxis.label.set_size(8)
c.grid(False)
d.grid(False)
family_header_c=c.inset_axes(CONFIG["heatmap"]["header_bounds"])
family_header_d=d.inset_axes(CONFIG["heatmap"]["header_bounds"])
family_header_palette=ListedColormap([NORMAL,JOINT_COLOR,DV_COLOR])
sns.heatmap(np.array([[0,1,2]]),ax=family_header_c,cmap=family_header_palette,vmin=0,vmax=2,
    cbar=False,xticklabels=False,yticklabels=False,linewidths=2,linecolor="white")
sns.heatmap(np.array([[0,1,2]]),ax=family_header_d,cmap=family_header_palette,vmin=0,vmax=2,
    cbar=False,xticklabels=False,yticklabels=False,linewidths=2,linecolor="white")
fig.savefig(F / "Figure_2.pdf")
fig.savefig(F / "Figure_2.svg")
fig.savefig(F / "Figure_2.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["4"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure 4: paired curves at left and method comparison at right.
# Effect scatter and marginal distributions follow QuantUMS Fig. 1b.
q = pd.read_csv(R / f"sgcb_{CONFIG['mean_record_id']}_seqc_qpcr_metrics.csv")
q = q.loc[q.reference_estimand.eq("geometric_mean")]
q = q.loc[pd.to_numeric(q.qpcr_abs_log2FC_threshold, errors="coerce").eq(1)].sort_values("AUC", ascending=False)
scores = pd.read_csv(R / f"sgcb_{CONFIG['mean_record_id']}_seqc_qpcr_gene_scores.csv")
own = scores.loc[scores.method.eq("SGCB")]
fig = plt.figure(figsize=(7.25, 5.0), layout="constrained")
grid = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1.35, 1])
scgrid = grid[0, 0].subgridspec(2, 2, width_ratios=[6, 1], height_ratios=[1, 5], wspace=.015, hspace=.015)
ax = fig.add_subplot(scgrid[1, 0])
ax.grid(True)
sns.scatterplot(data=own,x="log2FC_qpcr_B_over_A",y="log2FC",s=18,color=MEAN_COLOR,alpha=.74,edgecolor="white",linewidth=.25,ax=ax,legend=False)
ax.axline((0, 0), slope=1, lw=.85, ls="--", color=NORMAL)
ax.set(xlabel="TaqMan log₂(B/A)", ylabel="SGCB log₂ mean ratio")
ax.text(.04, .96, f"{len(own):,} genes\nSpearman ρ = {own.log2FC_qpcr_B_over_A.corr(own.log2FC, method='spearman'):.3f}", transform=ax.transAxes, va="top", fontsize=8)
top = fig.add_subplot(scgrid[0, 0], sharex=ax)
side = fig.add_subplot(scgrid[1, 1], sharey=ax)
sns.histplot(x=own.log2FC_qpcr_B_over_A,bins="fd",color=FIG_STYLE["soft"],element="step",fill=True,alpha=.9,linewidth=.7,ax=top)
sns.histplot(y=own.log2FC,bins="fd",color=FIG_STYLE["soft"],element="step",fill=True,alpha=.9,linewidth=.7,ax=side)
top.set_axis_off()
side.set_axis_off()
top.set_title("$\\mathbf{a}$")
ax = fig.add_subplot(grid[1, :])
ax.grid(axis="x")
auc_colors=q.method.map(FIG_STYLE["methods"]).to_numpy()
ax.hlines(np.arange(len(q)),q.AUC_low95,q.AUC_high95,colors=auc_colors,lw=CONFIG["detail"]["interval_width"])
ax.scatter(q.AUC,np.arange(len(q)),c=auc_colors,s=30,edgecolors="white",linewidths=.5,zorder=3)
ax.scatter(q.AUC_low95,np.arange(len(q)),c=auc_colors,marker="|",s=22,linewidths=1)
ax.scatter(q.AUC_high95,np.arange(len(q)),c=auc_colors,marker="|",s=22,linewidths=1)
position = np.flatnonzero(q.method.eq("SGCB"))[0]
ax.axhspan(position-.34,position+.34,color=MEAN_COLOR,alpha=.045,zorder=0)
sg = q.loc[q.method.eq("SGCB")].iloc[0]
ax.errorbar(sg.AUC, position, xerr=np.array([[sg.AUC - sg.AUC_low95], [sg.AUC_high95 - sg.AUC]]),
            fmt="D", color=MEAN_COLOR, capsize=CONFIG["detail"]["interval_caps"],
            elinewidth=CONFIG["detail"]["interval_width"],ms=5,mec="white",mew=.5,zorder=4)
auc_labels = np.frompyfunc(ax.text, 3, 1)(np.ones(len(q)), np.arange(len(q)), np.char.mod("%.3f", q.AUC))
plt.setp(auc_labels.tolist(), fontsize=8, ha="left", va="center", clip_on=False,
    transform=ax.get_yaxis_transform() + ScaledTranslation(
        plt.rcParams["figure.constrained_layout.w_pad"], 0, fig.dpi_scale_trans))
ax.set(yticks=np.arange(len(q)), yticklabels=q.method, xlim=(.68, 1), ylim=(len(q) - .45, -.9),
       xlabel="ROC AUC (95% CI)", title="$\\mathbf{c}$", xticks=[.7, .8, .9])
ax.get_yticklabels()[position].set_color(MEAN_COLOR)
ax.get_yticklabels()[position].set_fontweight("bold")
ax.text(.04, .98, "|qPCR log₂ fold change| ≥ 1", fontsize=7.0, va="top", transform=ax.transAxes)
ax = fig.add_subplot(grid[0, 1])
ax.grid(True)
pr = scores.assign(truth=scores.log2FC_qpcr_B_over_A.abs().ge(1), score=-np.log10(scores.pvalue.clip(lower=1e-300)))
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("SGCB"), "truth"], pr.loc[pr.method.eq("SGCB"), "score"])
ax.plot(rec, p, color=MEAN_COLOR, lw=1.7, label="SGCB", zorder=5)
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("DESeq2"), "truth"], pr.loc[pr.method.eq("DESeq2"), "score"])
ax.plot(rec, p, color=FIG_STYLE["methods"]["DESeq2"], lw=1.1, label="DESeq2")
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("edgeR"), "truth"], pr.loc[pr.method.eq("edgeR"), "score"])
ax.plot(rec, p, color=FIG_STYLE["methods"]["edgeR"], lw=1.1, label="edgeR")
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("limma"), "truth"], pr.loc[pr.method.eq("limma"), "score"])
ax.plot(rec, p, color=FIG_STYLE["methods"]["limma"], ls="--", lw=1.1, label="limma")
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("glmGamPoi"), "truth"], pr.loc[pr.method.eq("glmGamPoi"), "score"])
ax.plot(rec, p, color=FIG_STYLE["methods"]["glmGamPoi"], ls=":", lw=1.2, label="glmGamPoi")
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("NOISeq"), "truth"], pr.loc[pr.method.eq("NOISeq"), "score"])
ax.plot(rec, p, color=FIG_STYLE["methods"]["NOISeq"], lw=1.1, label="NOISeq")
p, rec, _ = precision_recall_curve(pr.loc[pr.method.eq("EBSeq"), "truth"], pr.loc[pr.method.eq("EBSeq"), "score"])
ax.plot(rec, p, color=FIG_STYLE["methods"]["EBSeq"], lw=1.1, label="EBSeq")
ax.axhline(pr.loc[pr.method.eq("SGCB"), "truth"].mean(), color=GREY, ls="--", lw=.7)
ax.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1.015), ylim=(.55, 1.015), title="$\\mathbf{b}$")
comparison_handles,comparison_labels=ax.get_legend_handles_labels()
fig.legend(comparison_handles,comparison_labels,frameon=False,fontsize=CONFIG["detail"]["legend_size"],
    ncols=4,loc="outside lower center",handlelength=1.7,columnspacing=1,labelspacing=.5)
fig.savefig(F / "Figure_4.pdf")
fig.savefig(F / "Figure_4.svg")
fig.savefig(F / "Figure_4.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["5"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure 5: cohort, functional context and complete gene distributions.
kirc = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_genes.csv")
stability = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_distribution_stability.csv")
examples = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_distribution_examples.csv")
features = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_example_features.csv").set_index("gene_id")
cov = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_covariate_genes.csv")
stability_summary = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_distribution_stability_summary.csv").set_index("channel")
dv_stability = stability_summary.loc["DV_only"]
enrich = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_enrichment.csv")
esum = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_enrichment_summary.csv")
counts = kirc.channel.value_counts().reindex(CHANNELS)
fig = plt.figure(figsize=(7.25, 7.1), layout="constrained")
outer = fig.add_gridspec(3,2,width_ratios=[1,1.38],height_ratios=[.55,1.10,1])
left=[outer[0,0],outer[1,0],outer[2,0]]
right=[outer[0,1],outer[1,1],outer[2,1]]
ax = fig.add_subplot(left[0])
ax.set_title("$\\mathbf{a}$")
ax.set_axis_off()
mosaic_counts={("DE","DV"):int(counts["DE_and_DV"]),("DE","No DV"):int(counts["DE_only"]),
    ("No DE","DV"):int(counts["DV_only"]),("No DE","No DV"):int(counts["neither"])}
mosaic_labels={("DE","DV"):f"Both\n{counts['DE_and_DV']:,}",("DE","No DV"):f"DE only\n{counts['DE_only']:,}",
    ("No DE","DV"):f"DV only\n{counts['DV_only']:,}",("No DE","No DV"):f"Neither\n{counts['neither']:,}"}
mosaic_colors={("DE","DV"):{"facecolor":JOINT_COLOR,"edgecolor":"white","linewidth":1.2},
    ("DE","No DV"):{"facecolor":FIG_STYLE["cohorts"]["normal"],"edgecolor":"white","linewidth":1.2},
    ("No DE","DV"):{"facecolor":FIG_STYLE["cohorts"]["tumor"],"edgecolor":"white","linewidth":1.2},
    ("No DE","No DV"):{"facecolor":FIG_STYLE["colors"]["other"],"edgecolor":"white","linewidth":1.2}}
mosaic(mosaic_counts,ax=ax,gap=.018,properties=mosaic_colors,labelizer=lambda key:mosaic_labels[key],axes_label=False,title=r"$\mathbf{a}$")
plt.setp(ax.texts,fontsize=8,color=INK)
ax.text(0,.98,"KIRC · 59 normal · 80 tumors",transform=ax.transAxes,va="top",fontsize=7.5)
ax.set(xlim=(0,1),ylim=(0,1.25))
ax = fig.add_subplot(left[1])
ax.grid(True)
sv = stability.loc[stability.channel.eq("DV_only") & stability.empirical_cv2_discovery.gt(0) & stability.empirical_cv2_heldout.gt(0)]
hm = ax.hexbin(np.log10(sv.empirical_cv2_discovery), np.log10(sv.empirical_cv2_heldout),
               cmap=sns.blend_palette(FIG_STYLE["density_scale"],as_cmap=True), norm=LogNorm(), **CONFIG["hexbin"])
ax.axline((0, 0), slope=1, color=NORMAL, ls="--", lw=.8)
ax.set(title="$\\mathbf{c}$", xlabel="Discovery tumor log₁₀ CV²", ylabel="Held-out tumor log₁₀ CV²")
ax.text(.03, .96, f"{len(sv):,} DV-only genes\n{int(dv_stability.discovery_tumors)} vs {int(dv_stability.heldout_tumors)} tumors\nρ = {dv_stability.rho_tumor_cv2:.3f}", va="top", transform=ax.transAxes, fontsize=7.5)
ax = fig.add_subplot(left[2])
ax.grid(True)
cvsel = cov.loc[cov.channel.eq("DV_only")].dropna(subset=["adjusted_log2_cv2_ratio", "unadjusted_log2_cv2_ratio"])
ax.scatter(cvsel.unadjusted_log2_cv2_ratio, cvsel.adjusted_log2_cv2_ratio, s=7, color=DV_COLOR, alpha=.4, linewidths=0)
ax.axline((0, 0), slope=1, color=NORMAL, ls="--", lw=.8)
ax.set(title="$\\mathbf{e}$", xlabel="Unadjusted moment log₂ CV² ratio", ylabel="Adjusted moment log₂ CV² ratio")
ax.text(.03, .97, f"{len(cvsel):,} DV-only genes\nAge, sex, immune proxy, purity", va="top", transform=ax.transAxes, fontsize=7.2)
zoom_limits=np.quantile(cvsel[["unadjusted_log2_cv2_ratio","adjusted_log2_cv2_ratio"]].to_numpy(),
    CONFIG["kirc_zoom"]["quantiles"])
zoom=ax.inset_axes(CONFIG["kirc_zoom"]["bounds"])
zoom.scatter(cvsel.unadjusted_log2_cv2_ratio,cvsel.adjusted_log2_cv2_ratio,
    s=3,color=DV_COLOR,alpha=.35,linewidths=0,rasterized=True)
zoom.axline((0,0),slope=1,color=MEAN_COLOR,ls="--",lw=.6)
zoom.set(xlim=zoom_limits,ylim=zoom_limits)
zoom.tick_params(labelsize=6.3,length=2)

ax = fig.add_subplot(right[0])
dv_enrich = enrich.loc[enrich.category.eq("DV_only") & enrich.qvalue.lt(CONFIG["alpha"])]
# Functional context and gene distributions follow Dreamlet Fig. 4.
selected = dv_enrich.set_index("term").reindex(CONFIG["kirc_enrichment_terms"]).reset_index()
term_positions = np.arange(len(selected))
term_axis_max = selected.fold_enrichment.max() * 1.35
points = ax.scatter(selected.fold_enrichment, term_positions, s=selected.k * .9 + 10,
                    c=-np.log10(selected.qvalue), cmap=sns.blend_palette(FIG_STYLE["density_scale"],as_cmap=True),
                    vmin=-np.log10(CONFIG["alpha"]), edgecolors=DV_COLOR, linewidths=.9)
term_labels=selected.term.map(CONFIG["kirc_enrichment_labels"])
ax.set(title="$\\mathbf{b}$", yticks=term_positions, yticklabels=term_labels,
       xlabel="Fold enrichment · DV-only set",
       xlim=(1, term_axis_max), ylim=(len(selected) - .4, -.6))
ax.axvline(1, color=GREY, ls=":", lw=.6)
np.frompyfunc(ax.text, 3, 1)(np.ones(len(selected)), term_positions,
    np.char.mod("q %.2g", selected.qvalue))
plt.setp(ax.texts,ha="right",va="center",fontsize=7.5,color=INK)
plt.setp(ax.texts,transform=ax.get_yaxis_transform()+ScaledTranslation(-3/72,0,fig.dpi_scale_trans))
ax.tick_params(axis="y",labelsize=8)
ax.tick_params(axis="y",length=0)
ax.spines["left"].set_visible(False)
ax.grid(axis="y")
gene_count_labels=np.frompyfunc(ax.text,3,1)(
    selected.fold_enrichment+(term_axis_max-1)*.035,term_positions,
    np.char.mod("%d genes",selected.k.to_numpy()))
plt.setp(gene_count_labels.tolist(),fontsize=7.0,va="center",ha="left",color=INK)
ax.get_yticklabels()[0].set_color(DV_COLOR)
ax.get_yticklabels()[0].set_fontweight("bold")
dv_terms = esum.loc[esum.category.eq("DV_only")].set_index("collection").significant_terms

genegrid = right[1].subgridspec(1, 3,width_ratios=[1.4,1,1])
ax1, ax2, ax3 = fig.add_subplot(genegrid[0]), fig.add_subplot(genegrid[1]), fig.add_subplot(genegrid[2])
ga = examples.loc[examples.gene_id.eq("SF3B1")]
gb = examples.loc[examples.gene_id.eq("FTL")]
gc = examples.loc[examples.gene_id.eq("GPX3")]
ga=ga.assign(log_abundance=np.log2(ga.abundance+.5))
gb=gb.assign(log_abundance=np.log2(gb.abundance+.5))
gc=gc.assign(log_abundance=np.log2(gc.abundance+.5))
cohort_order=["normal","tumor","heldout_tumor"]
sns.boxenplot(data=ga,x="log_abundance",y="group",hue="group",order=cohort_order,hue_order=cohort_order,
    palette=FIG_STYLE["cohorts"],legend=False,k_depth=3,showfliers=False,saturation=1,
    width=.62,linewidth=.85,linecolor=CONFIG["detail"]["outline"],ax=ax1)
sns.boxenplot(data=gb,x="log_abundance",y="group",hue="group",order=cohort_order,hue_order=cohort_order,
    palette=FIG_STYLE["cohorts"],legend=False,k_depth=3,showfliers=False,saturation=1,
    width=.62,linewidth=.85,linecolor=CONFIG["detail"]["outline"],ax=ax2)
sns.boxenplot(data=gc,x="log_abundance",y="group",hue="group",order=cohort_order,hue_order=cohort_order,
    palette=FIG_STYLE["cohorts"],legend=False,k_depth=3,showfliers=False,saturation=1,
    width=.62,linewidth=.85,linecolor=CONFIG["detail"]["outline"],ax=ax3)
sns.swarmplot(data=ga,x="log_abundance",y="group",order=cohort_order,color="#526F75",size=CONFIG["distribution"]["point_size"],alpha=.48,ax=ax1)
sns.swarmplot(data=gb,x="log_abundance",y="group",order=cohort_order,color="#526F75",size=CONFIG["distribution"]["point_size"],alpha=.48,ax=ax2)
sns.swarmplot(data=gc,x="log_abundance",y="group",order=cohort_order,color="#526F75",size=CONFIG["distribution"]["point_size"],alpha=.48,ax=ax3)
ax1.set(title=r"$\mathbf{d}$",xlabel="Log₂(CPM + 0.5)",ylabel="")
ax1.set_title(r"$\mathbf{d}$",loc="left",color=INK,fontweight="bold")
ax2.set(title="",xlabel="Log₂(CPM + 0.5)",ylabel="")
ax3.set(title="",xlabel="Log₂(CPM + 0.5)",ylabel="")
ax1.text(.97,.97,"SF3B1",transform=ax1.transAxes,ha="right",va="top",
    color=DV_COLOR,fontsize=8,fontweight="bold")
ax2.text(.97,.97,"FTL",transform=ax2.transAxes,ha="right",va="top",
    color=INK,fontsize=8,fontweight="bold")
ax3.text(.97,.97,"GPX3",transform=ax3.transAxes,ha="right",va="top",
    color=INK,fontsize=8,fontweight="bold")
plt.setp([ax1,ax2,ax3],ylim=(2.55,-.9),yticks=[0,1,2],yticklabels=[])
ax1.grid(axis="x")
ax2.grid(axis="x")
ax3.grid(axis="x")
ax1.text(.03,.96,f"DV q {features.loc['SF3B1','dv_padj']:.3g}",transform=ax1.transAxes,
    color=INK,fontsize=7.5,va="top")
ax1.text(.03,.86,f"Mean q {features.loc['SF3B1','padj']:.3g}",transform=ax1.transAxes,
    color=INK,fontsize=7.5,va="top")
ax2.text(.03, .96, CHANNEL_LABELS[features.loc["FTL", "channel"]], transform=ax2.transAxes,
    color=INK, fontsize=7.5, va="top")
ax3.text(.03, .96, CHANNEL_LABELS[features.loc["GPX3", "channel"]], transform=ax3.transAxes,
    color=INK, fontsize=7.5, va="top")
ax = fig.add_subplot(right[2])
cx = cov.set_index("gene_id").loc[["SF3B1", "FTL", "GPX3"]]
cov_summary = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_covariate_summary.csv").iloc[0]
cov_critical = t.ppf(.975, cov_summary.n_donors - 1)
ax.errorbar(cx.unadjusted_log2_cv2_ratio, np.arange(3) - .13, xerr=cov_critical * cx.unadjusted_se_log2_cv2,
            fmt="o", color=MEAN_COLOR,capsize=CONFIG["detail"]["interval_caps"],
            elinewidth=CONFIG["detail"]["interval_width"],ms=4.5,mec="white",mew=.5,label="Group-only")
ax.errorbar(cx.adjusted_log2_cv2_ratio, np.arange(3) + .13, xerr=cov_critical * cx.adjusted_se_log2_cv2,
            fmt="s", color=DV_COLOR,capsize=CONFIG["detail"]["interval_caps"],
            elinewidth=CONFIG["detail"]["interval_width"],ms=5,mec="white",mew=.5,label="Adjusted")
ax.axvline(0, color=NORMAL, ls=":", lw=.8)
ax.axhspan(-.38, .4, color=DV_COLOR, alpha=.055, zorder=0)
ax.set(title="$\\mathbf{f}$", yticks=np.arange(3), yticklabels=cx.index,
       xlabel="Moment log₂ CV² ratio (95% CI)", ylim=(2.8, -.75), xlim=(-3.5, 6.2), xticks=[-2, 0, 2, 4])
ax.get_yticklabels()[0].set_color(DV_COLOR)
ax.get_yticklabels()[0].set_fontweight("bold")
np.frompyfunc(ax.text, 3, 1)(np.repeat(6.1, 3), np.arange(3) + .13,
    np.char.mod("q = %.3g", cx.adjusted_qvalue.to_numpy()))
plt.setp(ax.texts, fontsize=7.8, va="center", ha="right", color=INK)
ax.legend(frameon=False, loc="lower right", fontsize=7.5, handlelength=1.3, ncols=2)
fig.legend(handles=[Patch(facecolor=FIG_STYLE["cohorts"]["normal"],edgecolor="#65868C",label="Normal (59)"),Patch(facecolor=FIG_STYLE["cohorts"]["tumor"],edgecolor="#65868C",label="Discovery tumor (80)"),Patch(facecolor=FIG_STYLE["cohorts"]["heldout_tumor"],edgecolor="#65868C",label="Held-out tumor (394)")],frameon=False,fontsize=7.5,loc="outside lower center",ncols=3)
fig.savefig(F / "Figure_5.pdf")
fig.savefig(F / "Figure_5.svg")
fig.savefig(F / "Figure_5.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["2"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S6 compares the same CV2 target, then shows why distribution choice matters.
sim = pd.read_csv(R / "current_simulation_metrics.csv")
sc = pd.read_csv(R / "current_simulation_scenarios.csv")
design_order = CONFIG["designs"]
emp = pd.concat([pd.read_csv(R / "moment_cv_abundance_n10_10_summary.csv").assign(design="10:10"),
    pd.read_csv(R / "moment_cv_abundance_n10_50_summary.csv").assign(design="10:50"),
    pd.read_csv(R / "moment_cv_abundance_n50_50_summary.csv").assign(design="50:50")])
ln = pd.concat([pd.read_csv(R / "simulation_lognormal_n10_10_summary.csv").assign(design="10:10"),
    pd.read_csv(R / "simulation_lognormal_n10_50_summary.csv").assign(design="10:50"),
    pd.read_csv(R / "simulation_lognormal_n50_50_summary.csv").assign(design="50:50")])
ln["sensitivity"] = ln.TP / 1200
current_cv = api_metrics.loc[api_metrics.panel.eq("channel") & api_metrics.endpoint.eq("DV")].copy()
current_cv["sensitivity"] = current_cv.power
cv_comparison = pd.concat([current_cv.assign(method="SGCB"),
    sim.loc[sim.target.eq("CV2")].assign(method="Gene-wise GG"),
    emp.assign(method=emp.test.map({"p_delta": "Empirical delta", "p_jackknife": "Jackknife"})),
    ln.assign(method="Lognormal F")], ignore_index=True)
cv_comparison[["design", "method", "TP", "FP", "FDP", "sensitivity"]].to_csv(
    R / "publication_CV2_comparison.csv", index=False, encoding="utf-8")
methods = ["SGCB", "Gene-wise GG", "Empirical delta", "Jackknife", "Lognormal F"]
method_colors = [DV_COLOR, GREY, NORMAL, MEAN_COLOR, JOINT_COLOR]
power = cv_comparison.pivot(index="design", columns="method", values="sensitivity").reindex(design_order)[methods]
fdp = cv_comparison.pivot(index="design", columns="method", values="FDP").reindex(design_order)[methods]
xp = (np.arange(3)[:, None] + np.linspace(-.32,.32,len(methods))[None, :]).ravel()
colors = np.tile(method_colors, 3)
fig = plt.figure(figsize=(7.8, 7.05), layout="constrained")
grid = fig.add_gridspec(2, 2, height_ratios=[.85, 1.15])
ax = fig.add_subplot(grid[0, 0])
ax.bar(xp, power.to_numpy().ravel(), width=.15, color=colors,edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
ax.set(title="$\\mathbf{a}$  Sensitivity for CV²", xticks=np.arange(3), xticklabels=design_order,
       xlabel="Samples per group", ylabel="True positives / 1,200", ylim=(0, 1.05))
ax = fig.add_subplot(grid[0, 1])
ax.bar(xp, fdp.to_numpy().ravel(), width=.15, color=colors,edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
ax.axhline(.05, color=RED, ls=":", lw=1)
ax.text(2.45, .055, "5%", ha="right", va="bottom", color=RED, fontsize=8)
ax.set(title="$\\mathbf{b}$  False discoveries", xticks=np.arange(3), xticklabels=design_order,
       xlabel="Samples per group", ylabel="False / all discoveries", ylim=(0, .37))
ax = fig.add_subplot(grid[1, 0])
rate = api_scenarios.loc[api_scenarios.panel.eq("channel") & api_scenarios.design.eq("50:50")].pivot(
    index="scenario", columns="endpoint", values="calls").reindex(SCENARIOS)[["DE", "DV"]] / 600
sns.heatmap(rate,ax=ax,vmin=0,vmax=1,
    cmap=sns.blend_palette(FIG_STYLE["error_scale"],as_cmap=True),
    annot=np.char.mod("%.1f",rate.to_numpy()*100),fmt="",annot_kws={"fontsize":CONFIG["heatmap"]["font_size"]},
    linewidths=CONFIG["heatmap"]["cell_gap"],linecolor="white",cbar=False,
    xticklabels=["DE: arithmetic\nmean","DV: CV²"],yticklabels=SCENARIO_LABELS)
ax.set(title=r"$\mathbf{c}$  Mean and CV² responses (%)",xlabel="",ylabel="")
ax.tick_params(length=0,labelrotation=0)
ax.spines[["left", "bottom"]].set_visible(False)
ax.grid(False)
response_header=ax.inset_axes(CONFIG["heatmap"]["header_bounds"])
sns.heatmap(np.array([[0,1]]),ax=response_header,cmap=ListedColormap([MEAN_COLOR,DV_COLOR]),
    vmin=0,vmax=1,cbar=False,xticklabels=False,yticklabels=False,linewidths=2,linecolor="white")
ax = fig.add_subplot(grid[1, 1])
mom = pd.read_csv(R / "CV2_estimator_comparison.csv")
newmom = pd.read_csv(R / "current_CV2_estimation.csv")
rmse_order = ["null_gg", "null_lognormal", "null_contamination", "null_mixture"]
rmse = mom.pivot(index="scenario", columns="estimator", values="log_CV2_RMSE").reindex(rmse_order)
newrmse = newmom.loc[newmom.design.eq("50:50")].set_index("scenario").reindex(rmse_order)
ax.scatter(newrmse.log_CV2_RMSE, np.arange(4) - .16, color=GREY,marker="D",s=34,edgecolors=CONFIG["detail"]["outline"],linewidths=.5,label="Gene-wise GG")
ax.scatter(rmse.sample_moments, np.arange(4), color=NORMAL, marker="s", s=30,edgecolors="white",linewidths=.5,label="Empirical moments")
ax.scatter(rmse.lognormal, np.arange(4) + .16, color=JOINT_COLOR, marker="^", s=36,edgecolors="white",linewidths=.5,label="Lognormal fit")
ax.set(title="$\\mathbf{d}$  Gene-wise CV² estimation", yticks=np.arange(4),
       yticklabels=["Gamma", "Lognormal", "Contaminated", "Mixture"], xlabel="RMSE of log CV²", xlim=(0, 2.03), ylim=(3.5, -.55))
ax.legend(frameon=False, loc="lower right", fontsize=7.5)
fig.legend(handles=[Line2D([], [], marker="s", ls="", color=DV_COLOR, label="SGCB"),
    Line2D([], [], marker="s", ls="", color=GREY, label="Gene-wise GG"),
    Line2D([], [], marker="s", ls="", color=NORMAL, label="Empirical delta"),
    Line2D([], [], marker="s", ls="", color=MEAN_COLOR, label="Jackknife"),
    Line2D([], [], marker="s", ls="", color=JOINT_COLOR, label="Lognormal F")],
    loc="outside lower center", ncols=3, frameon=False, fontsize=8)
fig.savefig(F / "Figure_S6.pdf")
fig.savefig(F / "Figure_S6.svg")
fig.savefig(F / "Figure_S6.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["1"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S1 separates zero representation, regularization and covariance estimation.
factorial = pd.read_csv(R / "zero_shape_factorial_summary.csv", keep_default_na=False)
factorial["method"] = factorial.zero_model + " / " + factorial.shape_prior
factorial_order = ["Offset continuous / None", "Offset continuous / N(0,1)", "Zero mass / None", "Zero mass / N(0,1)"]
zero_order = ["None", "(0,20%]", "(20,50%]", "(50,100%]"]
conv = factorial.pivot(index="zero_class", columns="method", values="converged_both").reindex(zero_order)[factorial_order]
usable = factorial.pivot(index="zero_class", columns="method", values="eligible_finite_cv").reindex(zero_order)[factorial_order]
fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.7), layout="constrained")
xp = (np.arange(4)[:, None] + np.array([-.27, -.09, .09, .27])[None, :]).ravel()
component_colors = [GREY,NORMAL,CONFIG["figure_palettes"]["1"]["soft"],DV_COLOR]
axes[0, 0].bar(xp, conv.to_numpy().ravel(), width=.15, color=np.tile(component_colors, 4),edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
axes[0, 0].set(title="$\\mathbf{a}$  Fitting convergence", xticks=np.arange(4), xticklabels=["0", "≤20%", "20–50%", ">50%"],
               xlabel="Fraction of zero observations", ylabel="Both groups converged", ylim=(0, 1.06))
axes[0, 1].bar(xp, usable.to_numpy().ravel(), width=.15, color=np.tile(component_colors, 4),edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
axes[0,0].grid(axis="y")
axes[0,1].grid(axis="y")
axes[0, 1].set(title="$\\mathbf{b}$  Finite CV² inference", xticks=np.arange(4), xticklabels=["0", "≤20%", "20–50%", ">50%"],
               xlabel="Fraction of zero observations", ylabel="Eligible fraction", ylim=(0, 1.06))
diag = pd.concat([pd.read_csv(R / "penalized_covariance_abundance_n10_10_summary.csv").assign(design="10:10"),
    pd.read_csv(R / "penalized_covariance_abundance_n10_50_summary.csv").assign(design="10:50"),
    pd.read_csv(R / "penalized_covariance_abundance_n50_50_summary.csv").assign(design="50:50")])
old = diag.loc[diag.covariance.eq("Posterior Hessian")].set_index("design").reindex(design_order)
new = diag.loc[diag.covariance.eq("Expected information")].set_index("design").reindex(design_order)
axes[1, 0].hlines(np.arange(3), old.sensitivity, new.sensitivity, color=GREY, lw=1.5)
axes[1, 0].scatter(old.sensitivity, np.arange(3), color=NORMAL, marker="s",edgecolors="white",linewidths=.5,label="Inverse penalized Hessian",s=34)
axes[1, 0].scatter(new.sensitivity, np.arange(3), color=DV_COLOR,edgecolors="white",linewidths=.5,label="Expected sampling covariance",s=42)
axes[1, 0].set(title="$\\mathbf{c}$  Sampling uncertainty", yticks=np.arange(3), yticklabels=design_order,
               xlabel="CV² sensitivity", ylabel="Samples per group", xlim=(-.03, 1.01), ylim=(-.5, 3.15))
axes[1, 0].legend(frameon=False, loc="upper left", fontsize=7.8)
null = pd.read_csv(R / "expected_information_null_summary.csv")
axes[1, 1].scatter(null.split, null.raw_p005, color=DV_COLOR, s=32)
axes[1, 1].axhline(.05, color=RED, ls=":", lw=1)
axes[1, 1].set(title="$\\mathbf{d}$  Gene-wise null diagnostic", xlabel="Partition", ylabel="Unadjusted p < 0.05 fraction",
               xticks=np.arange(1, 11), ylim=(0, .07))
axes[1, 1].text(.04, .94, "0 BH discoveries in each of 10 partitions", va="top", transform=axes[1, 1].transAxes, fontsize=8)
fig.legend(handles=[Line2D([], [], marker="s", ls="", color=component_colors[0], label="Offset, no prior"),
    Line2D([], [], marker="s", ls="", color=component_colors[1], label="Offset, Q prior"),
    Line2D([], [], marker="s", ls="", color=component_colors[2], label="Zero mass, no prior"),
    Line2D([], [], marker="s", ls="", color=component_colors[3], label="Zero mass, Q prior")],
    frameon=False, ncols=2, loc="outside lower center", fontsize=CONFIG["detail"]["legend_size"],columnspacing=1.8,labelspacing=.5)
fig.savefig(F / "Figure_S1.pdf")
fig.savefig(F / "Figure_S1.svg")
fig.savefig(F / "Figure_S1.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["5"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S2 uses current discovery memberships and the fixed score definitions.
NORMAL=FIG_STYLE["colors"]["both"]
fig = plt.figure(figsize=(7.8, 6.3), layout="constrained")
grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
ax = fig.add_subplot(grid[0, 0])
tcga = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_heldout_survival.csv").query("model == 'clinical'").set_index("term").loc[CHANNELS[:3]]
ext = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_external_survival.csv").query("model == 'clinical'").set_index("term").loc[CHANNELS[:3]]
yy = np.arange(3)
ax.errorbar(tcga.estimate, yy - .12, xerr=np.vstack([tcga.estimate - tcga["conf.low"], tcga["conf.high"] - tcga.estimate]),
            fmt="o", color=NORMAL,capsize=CONFIG["detail"]["interval_caps"],elinewidth=CONFIG["detail"]["interval_width"],mec="white",mew=.5,label="TCGA held-out (385)")
ax.errorbar(ext.estimate, yy + .12, xerr=np.vstack([ext.estimate - ext["conf.low"], ext["conf.high"] - ext.estimate]),
            fmt="s", color=TUMOR,capsize=CONFIG["detail"]["interval_caps"],elinewidth=CONFIG["detail"]["interval_width"],mec="white",mew=.5,label="External (99)")
ax.axvline(1, color=GREY, ls="--", lw=.8)
ax.set(yticks=yy, yticklabels=["DE only", "DV only", "Both"], xlabel="Adjusted HR per score SD (95% CI)",
       xlim=(min(tcga["conf.low"].min(), ext["conf.low"].min()) * .9,
             max(tcga["conf.high"].max(), ext["conf.high"].max()) * 1.05),
       ylim=(2.9, -.35), title="$\\mathbf{a}$  SGCB discovery gene sets")
cohort_handles,cohort_labels=ax.get_legend_handles_labels()
ax.grid(axis="x")
ax = fig.add_subplot(grid[0, 1])
alt = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_alternative_scores_survival.csv")
alte = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_external_alternative_scores_survival.csv")
ax.errorbar(alt.estimate, np.arange(2) - .1, xerr=np.vstack([alt.estimate - alt["conf.low"], alt["conf.high"] - alt.estimate]),
            fmt="o", color=NORMAL,capsize=CONFIG["detail"]["interval_caps"],elinewidth=CONFIG["detail"]["interval_width"],mec="white",mew=.5,label="TCGA held-out")
ax.errorbar(alte.estimate, np.arange(2) + .1, xerr=np.vstack([alte.estimate - alte["conf.low"], alte["conf.high"] - alte.estimate]),
            fmt="s", color=TUMOR,capsize=CONFIG["detail"]["interval_caps"],elinewidth=CONFIG["detail"]["interval_width"],mec="white",mew=.5,label="External")
ax.axvline(1, color=GREY, ls="--", lw=.8)
ax.set(yticks=[0, 1], yticklabels=["Direction-weighted\nmean", "Absolute\ndeviation"], ylim=(1.65, -.5),
       xlabel="Adjusted HR per score SD (95% CI)", title="$\\mathbf{b}$  Fixed score sensitivities")
ax.grid(axis="x")
ax = fig.add_subplot(grid[1, :])
km = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_heldout_KM.csv")
km_test = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_heldout_KM_test.csv").iloc[0]
lo, hi = km.loc[km.strata.eq("DV_high=Low")], km.loc[km.strata.eq("DV_high=High")]
ax.step(np.r_[0, lo.time / 365.25], np.r_[1, lo.estimate], where="post", color=NORMAL, label=f"Low ({int(km_test.n_low)})")
ax.step(np.r_[0, hi.time / 365.25], np.r_[1, hi.estimate], where="post", color=TUMOR, label=f"High ({int(km_test.n_high)})")
ax.fill_between(lo.time / 365.25, lo["conf.low"], lo["conf.high"], step="post", color=NORMAL, alpha=.16, linewidth=0)
ax.fill_between(hi.time / 365.25, hi["conf.low"], hi["conf.high"], step="post", color=TUMOR, alpha=.16, linewidth=0)
ax.set(xlabel="Years", ylabel="Overall survival", ylim=(0, 1.05), title="$\\mathbf{c}$  Discovery-derived DV-only cutoff")
ax.text(.035, .14, f"High vs Low HR {km_test.estimate:.3f}\nLog-rank p = {km_test.logrank_p:.3f}", transform=ax.transAxes, fontsize=9)
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y")
fig.legend(cohort_handles,cohort_labels,frameon=False,loc="outside upper center",ncols=2,fontsize=CONFIG["detail"]["legend_size"])
fig.savefig(F / "Figure_S2.pdf")
fig.savefig(F / "Figure_S2.svg")
fig.savefig(F / "Figure_S2.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["6"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S3 uses the three current disease applications and their fixed covariate models.
sle_models = pd.read_csv(R / f"sgcb_{RECORD_ID}_sle_classical_monocytes_2021_covariate_summary.csv").set_index("model")
liver_models = pd.read_csv(R / f"sgcb_{RECORD_ID}_masld_liver_2023_covariate_summary.csv").set_index("model")
kirc_model = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_covariate_summary.csv").iloc[0]
kirc_covariates = pd.read_csv(R / f"sgcb_{RECORD_ID}_kirc_covariate_genes.csv")
kirc_selected = kirc_covariates.loc[kirc_covariates.population_DV &
    np.isfinite(kirc_covariates.primary_GG_effect) &
    np.isfinite(kirc_covariates.adjusted_log2_cv2_ratio) &
    np.isfinite(kirc_covariates.unadjusted_log2_cv2_ratio)].copy()
kirc_selected["same_direction"] = np.sign(kirc_selected.adjusted_log2_cv2_ratio).eq(
    np.sign(kirc_selected.primary_GG_effect))
sle_model = sle_models.loc["clinical"]
liver_model = liver_models.loc["clinical_and_composition"]
adjustment = pd.DataFrame({
    "dataset": ["sle_classical_monocytes_2021", "kirc", "masld_liver_2023"],
    "covariates": ["age + sex", "clinical covariates", "fibrosis + composition proxies"],
    "unadjusted_moment_calls": [sle_models.loc["unadjusted","calls"],
        kirc_model.unadjusted_moment_calls, liver_models.loc["unadjusted","moment_calls"]],
    "adjusted_calls": [sle_model.calls, kirc_model.adjusted_calls, liver_model.moment_calls],
    "population_calls": [sle_model.SGCB_selected, kirc_model.population_calls, liver_model.population_selected],
    "selected_available": [sle_model.selected_available, len(kirc_selected), liver_model.selected_available],
    "selected_same_direction_significant": [sle_model.selected_same_direction_significant,
        (kirc_selected.adjusted_qvalue.lt(.05) & kirc_selected.same_direction).sum(),
        liver_model.selected_same_direction_significant],
    "selected_direction_retained": [sle_model.selected_sign_retained,
        kirc_selected.same_direction.sum(), liver_model.selected_direction_retained]})
adjustment["selected_sign_retained"] = adjustment.selected_direction_retained / adjustment.selected_available
adjustment["selected_adjusted_fraction"] = adjustment.selected_same_direction_significant / adjustment.selected_available
adjustment.to_csv(R / "publication_CV2_adjustment.csv", index=False, encoding="utf-8")
labels = ["SLE monocytes", "KIRC", "MASLD liver"]
fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.9), layout="constrained")
axes[0].barh(np.arange(3)-.16, adjustment.unadjusted_moment_calls/1000,
    height=.28, color=GREY,edgecolor=MEAN_COLOR,linewidth=.5,label="Group only")
axes[0].barh(np.arange(3)+.16, adjustment.adjusted_calls/1000,
    height=.28, color=DV_COLOR,edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"],label="Adjusted")
axes[0].set(title="$\\mathbf{a}$  Conditional CV² screen", yticks=np.arange(3), yticklabels=labels,
    xlabel="Discoveries (thousands)", ylim=(2.6,-.6))
axes[0].legend(frameon=False, fontsize=7.5, loc="lower right")
bars = axes[1].barh(np.arange(3), adjustment.selected_adjusted_fraction, height=.55, color=DV_COLOR,edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
axes[1].bar_label(bars, labels=(adjustment.selected_same_direction_significant.astype(int).astype(str)+
    " / "+adjustment.selected_available.astype(int).astype(str)).to_numpy(), padding=3, fontsize=7.4)
axes[1].set(title="$\\mathbf{b}$  SGCB DV discoveries", yticks=np.arange(3), yticklabels=[],
    xlabel="Same direction\nAdjusted q < 0.05", xlim=(0,1.25),xticks=[0,.5,1],ylim=(2.6,-.6))
axes[2].scatter(adjustment.selected_sign_retained, np.arange(3), color=DV_COLOR, s=38,edgecolors="white",linewidths=.5)
axes[2].axvline(1, color=GREY, ls=":", lw=.7)
axes[2].set(title="$\\mathbf{c}$  Effect direction", yticks=np.arange(3), yticklabels=[],
    xlabel="Direction retained\n(fraction)", xlim=(0,1.03), ylim=(2.6,-.6))
axes[0].grid(axis="x")
axes[1].grid(axis="x")
axes[2].grid(axis="x")
fig.savefig(F / "Figure_S3.pdf")
fig.savefig(F / "Figure_S3.svg")
fig.savefig(F / "Figure_S3.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["2"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S4 compares dedicated tests against the separate log-spread truth.
dedicated = pd.concat([
    pd.read_csv(R / "simulation_diffVar_n10_10_summary.csv"),
    pd.read_csv(R / "simulation_diffVar_n10_50_summary.csv"),
    pd.read_csv(R / "simulation_diffVar_n50_50_summary.csv"),
    pd.read_csv(R / "simulation_clrDV_n10_10_summary.csv"),
    pd.read_csv(R / "simulation_clrDV_n10_50_summary.csv"),
    pd.read_csv(R / "simulation_clrDV_n50_50_summary.csv"),
    pd.read_csv(R / "simulation_fligner_n10_10_summary.csv"),
    pd.read_csv(R / "simulation_fligner_n10_50_summary.csv"),
    pd.read_csv(R / "simulation_fligner_n50_50_summary.csv")], ignore_index=True)
dedicated["design"] = dedicated.n0.astype(str) + ":" + dedicated.n1.astype(str)
all_tests = pd.read_csv(R / "current_simulation_test_results.csv")
bf = all_tests.loc[all_tests.target.eq("Log spread")].copy()
bf["CV2_TP"] = bf.called & bf.dv_truth
bf["CV2_FP"] = bf.called & ~bf.dv_truth
bf_cv = bf.groupby("design", sort=False)[["CV2_TP", "CV2_FP"]].sum().reset_index()
bf_summary = sim.loc[sim.target.eq("Log spread")].merge(bf_cv, on="design").assign(method="Brown–Forsythe")
dedicated = pd.concat([dedicated, bf_summary], ignore_index=True)
dedicated["method"] = dedicated.method.replace({"diffVar": "DiffVar", "clrDV": "clrDV", "fligner": "Fligner–Killeen"})
spread_alternatives = dedicated.design.map(bf.groupby("design").truth.sum())
dedicated["FN"] = spread_alternatives - dedicated.TP
dedicated["TN"] = dedicated.genes - spread_alternatives - dedicated.FP
dedicated["sensitivity"] = dedicated.TP / spread_alternatives
dedicated["precision"] = dedicated.TP / (dedicated.TP + dedicated.FP).clip(lower=1)
dedicated["F1"] = 2 * dedicated.TP / (2 * dedicated.TP + dedicated.FP + dedicated.FN)
dedicated["specificity"] = dedicated.TN / (dedicated.genes - spread_alternatives)
native = pd.concat([
    pd.read_csv(R / "simulation_clrDV_n10_10_genes.csv"),
    pd.read_csv(R / "simulation_clrDV_n10_50_genes.csv"),
    pd.read_csv(R / "simulation_clrDV_n50_50_genes.csv")], ignore_index=True)
native["design"] = native.n0.astype(str) + ":" + native.n1.astype(str)
native["native_BY_TP_logspread"] = native.q_original_BY.lt(.05) & native.log_spread_truth
native["native_BY_FP_logspread"] = native.q_original_BY.lt(.05) & ~native.log_spread_truth
native["native_BY_TP_CV2"] = native.q_original_BY.lt(.05) & native.dv_truth
native["native_BY_FP_CV2"] = native.q_original_BY.lt(.05) & ~native.dv_truth
native_summary = native.groupby(["method", "design"], sort=False)[[
    "native_BY_TP_logspread", "native_BY_FP_logspread", "native_BY_TP_CV2", "native_BY_FP_CV2"]].sum().reset_index()
native_summary["native_BY_FDP_logspread"] = native_summary.native_BY_FP_logspread / (
    native_summary.native_BY_TP_logspread + native_summary.native_BY_FP_logspread).clip(lower=1)
dedicated = dedicated.merge(native_summary, on=["method", "design"], how="left")
dedicated[["n0", "n1"]] = dedicated.design.str.split(":", expand=True).astype(int)
dedicated.to_csv(R / "publication_DV_comparators.csv", index=False, encoding="utf-8")
equal_moments = pd.read_csv(R / "equal_mean_cv2_distinct_log_dispersion.csv")
grid_x = np.linspace(.015, 50, 800)
alpha, beta, gamma = equal_moments.alpha.to_numpy()[:, None], equal_moments.beta.to_numpy()[:, None], equal_moments.gamma.to_numpy()[:, None]
log_density = np.log(gamma) - gamma * alpha * np.log(beta) - gammaln(alpha)
log_density = log_density + (gamma * alpha - 1) * np.log(grid_x)[None, :] - (grid_x[None, :] / beta) ** gamma
fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.4), layout="constrained")
axes[0, 0].plot(grid_x, np.exp(log_density[0]), color=NORMAL,lw=1.7,label="GG 1")
axes[0, 0].plot(grid_x, np.exp(log_density[1]), color=TUMOR,lw=1.7,ls="--",label="GG 2")
axes[0,0].fill_between(grid_x,0,np.exp(log_density[0]),color=NORMAL,alpha=.07)
axes[0,0].fill_between(grid_x,0,np.exp(log_density[1]),color=TUMOR,alpha=.07)
axes[0, 0].set(title="$\\mathbf{a}$  Equal abundance moments", xlabel="Abundance", ylabel="Density", xlim=(0, 40))
axes[0, 0].text(.97, .94, "Mean = 10\nCV² = 0.5", ha="right", va="top", transform=axes[0, 0].transAxes)
axes[0, 0].legend(frameon=False, loc="upper center")
axes[0, 1].scatter([-.09, .91], equal_moments.iloc[0][["log_variance", "mean_abs_log_deviation_from_median"]].to_numpy(dtype=float), color=NORMAL, s=36)
axes[0, 1].scatter([.09, 1.09], equal_moments.iloc[1][["log_variance", "mean_abs_log_deviation_from_median"]].to_numpy(dtype=float), color=TUMOR, s=36)
axes[0, 1].set(title="$\\mathbf{b}$  Different log spread", ylabel="Population quantity", xticks=[0, 1],
               xticklabels=["Log variance", "Absolute log deviation"], xlim=(-.4, 1.4), ylim=(.45, 1.1))
spread_methods = ["Brown–Forsythe", "DiffVar", "clrDV", "Fligner–Killeen"]
spread_colors = [NORMAL, MEAN_COLOR, JOINT_COLOR, TUMOR]
sp = dedicated.pivot(index="design", columns="method", values="sensitivity").reindex(CONFIG["designs"])[spread_methods]
sf = dedicated.pivot(index="design", columns="method", values="FDP").reindex(CONFIG["designs"])[spread_methods]
xp = (np.arange(3)[:, None] + np.array([-.27, -.09, .09, .27])[None, :]).ravel()
axes[1, 0].bar(xp, sp.to_numpy().ravel(), color=np.tile(spread_colors, 3), width=.15,edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
native_plot = native_summary.set_index("design").reindex(CONFIG["designs"])
axes[1, 0].scatter(np.arange(3) + .09, native_plot.native_BY_TP_logspread / 2400,
                   marker="D", facecolors="white", edgecolors=JOINT_COLOR,linewidths=1.4,s=34,zorder=4)
axes[1, 0].set(title="$\\mathbf{c}$  Log-spread sensitivity", ylabel="TP / 2,400 alternatives", xlabel="Samples per group",
               xticks=np.arange(3), xticklabels=CONFIG["designs"], ylim=(-.03, 1.03))
axes[1, 1].scatter(xp, sf.to_numpy().ravel(), c=np.tile(spread_colors, 3), s=33)
axes[1, 1].scatter(np.arange(3) + .09, native_plot.native_BY_FDP_logspread,
                   marker="D", facecolors="white", edgecolors=JOINT_COLOR,linewidths=1.4,s=34,zorder=4)
axes[1, 1].axhline(.05, color=RED, ls=":", lw=.8)
axes[1, 1].set(title="$\\mathbf{d}$  False discoveries", ylabel="Realized FDP", xlabel="Samples per group",
               xticks=np.arange(3), xticklabels=CONFIG["designs"],
               ylim=(-.02, max(.15, sf.max().max() * 1.12, native_plot.native_BY_FDP_logspread.max() * 1.12)))
spread_handles = [Line2D([], [], color=NORMAL, lw=3, label="Brown–Forsythe"),
    Line2D([], [], color=MEAN_COLOR, lw=3, label="DiffVar"), Line2D([], [], color=JOINT_COLOR, lw=3, label="clrDV"),
    Line2D([], [], color=TUMOR, lw=3, label="Fligner–Killeen"),
    Line2D([], [], color=JOINT_COLOR, ls="none", marker="D", markerfacecolor="white", label="clrDV native BY")]
fig.legend(handles=spread_handles, frameon=False, ncols=3, loc="outside lower center", fontsize=8)
fig.savefig(F / "Figure_S4.pdf")
fig.savefig(F / "Figure_S4.svg")
fig.savefig(F / "Figure_S4.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["4"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S5 reports the final public API, fixed-panel scaling and actual memory.
resource_labels = pd.Index("sgcb_"+RECORD_ID+"_"+pd.Series(CONFIG["resource_datasets"]),name="result_label")
resources = pd.concat(map(pd.read_csv,str(R)+"/"+resource_labels+"_summary.csv"),
    keys=resource_labels,names=["result_label","record"]).reset_index().drop(columns="record")
resources = resources.rename(columns={"elapsed_seconds":"api_seconds","cpu_seconds":"api_cpu_seconds",
    "process_peak_working_set_bytes":"api_peak_working_set_bytes",
    "process_peak_pagefile_bytes":"api_peak_pagefile_bytes"})
process_paths = str(R)+"/"+resource_labels+"_process_runtime.csv"
process_resources = pd.concat(map(pd.read_csv,process_paths),ignore_index=True)
resources = resources.merge(process_resources[["stage","elapsed_seconds","peak_working_set_bytes",
    "threads_per_process","execution_outcome"]].rename(columns={"stage":"result_label",
    "elapsed_seconds":"process_seconds","peak_working_set_bytes":"process_peak_working_set_bytes"}),
    on="result_label",how="left",validate="one_to_one")
resources["dataset"] = resources.result_label.str.removeprefix("sgcb_"+RECORD_ID+"_")
resources["api_peak_GiB"] = resources.api_peak_working_set_bytes/1024**3
resources["process_peak_GiB"] = resources.process_peak_working_set_bytes/1024**3
resources["api_minutes"] = resources.api_seconds/60
resources["samples"] = resources.n_ctrl+resources.n_treat
resources[["result_label","dataset","version","n_genes","n_ctrl","n_treat","api_seconds",
    "api_cpu_seconds","api_peak_GiB","process_seconds",
    "threads_per_process","execution_outcome"]].to_csv(
        R / "publication_current_resources.csv",index=False,encoding="utf-8")
sim_resources = resources.loc[resources.dataset.str.startswith(("channel_n","zero_n"))].copy()
sim_resources["panel"] = sim_resources.dataset.str.extract(r"^(channel|zero)_",expand=False)
positive_resources = sim_resources.loc[sim_resources.panel.eq("channel")].sort_values("samples")
zero_resources = sim_resources.loc[sim_resources.panel.eq("zero")].sort_values("samples")
cohort_resources = resources.set_index("dataset").loc[CONFIG["resource_cohorts"]]
fig,axes = plt.subplots(2,2,figsize=(7.25,6.0),layout="constrained")
axes[0,0].plot(positive_resources.samples,positive_resources.api_seconds,"o-",color=MEAN_COLOR,label="Positive (5,400 genes)")
axes[0,0].plot(zero_resources.samples,zero_resources.api_seconds,"s-",color=DV_COLOR,label="Observed zeros (4,000 genes)")
axes[0,0].set(title="$\\mathbf{a}$  Fixed-panel sample scaling",xlabel="Samples (two groups)",ylabel="API wall time (s)",xticks=[20,60,100])
axes[0,0].legend(frameon=False,fontsize=7.5)
axes[0,1].plot(positive_resources.samples,positive_resources.api_peak_GiB,"o-",color=MEAN_COLOR)
axes[0,1].plot(zero_resources.samples,zero_resources.api_peak_GiB,"s-",color=DV_COLOR)
axes[0,1].set(title="$\\mathbf{b}$  Fixed-panel memory",xlabel="Samples (two groups)",ylabel="Peak working set (GiB)",xticks=[20,60,100])
bars=axes[1,0].bar(np.arange(len(cohort_resources)),cohort_resources.api_minutes,color=MEAN_COLOR,width=.64,edgecolor="white",linewidth=CONFIG["detail"]["bar_edge_width"])
axes[1,0].bar_label(bars,fmt="%.1f",padding=4,fontsize=8.2)
axes[1,0].set(title="$\\mathbf{c}$  Complete cohort calls",ylabel="API wall time (min)",
    xticks=np.arange(len(cohort_resources)),xticklabels=CONFIG["resource_cohort_labels"],
    ylim=(0,cohort_resources.api_minutes.max()*1.2))
axes[1,0].tick_params(axis="x",labelsize=7.5)
resource_points=axes[1,1].scatter(resources.n_genes,resources.api_peak_GiB,
    s=CONFIG["resource_marker"]["offset"]+CONFIG["resource_marker"]["per_sample"]*resources.samples,
    color=DV_COLOR,alpha=.76,edgecolors="white",linewidths=.6)
axes[1,1].set(title=r"$\mathbf{d}$  All "+str(len(resources))+" complete calls",xlabel="Genes tested",
    ylabel="Peak working set (GiB)")
size_handles,size_labels=resource_points.legend_elements(prop="sizes",num=3,
    func=lambda area:(area-CONFIG["resource_marker"]["offset"])/CONFIG["resource_marker"]["per_sample"],
    color=DV_COLOR,markeredgecolor="white",markeredgewidth=.6,alpha=.76)
axes[1,1].legend(size_handles,size_labels,title="Samples",frameon=False,fontsize=7,title_fontsize=7.5,loc="upper right")
axes[0,0].grid(axis="y")
axes[0,1].grid(axis="y")
axes[1,0].grid(axis="y")
fig.savefig(F / "Figure_S5.pdf")
fig.savefig(F / "Figure_S5.svg")
fig.savefig(F / "Figure_S5.png")
plt.close(fig)


FIG_STYLE=CONFIG["figure_palettes"]["3"]
MEAN_COLOR,DV_COLOR,JOINT_COLOR=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
NORMAL,TUMOR,GREY=FIG_STYLE["colors"]["normal"],FIG_STYLE["colors"]["tumor"],FIG_STYLE["colors"]["reference"]
CHANNEL_COLORS={"DE_only":MEAN_COLOR,"DV_only":DV_COLOR,"DE_and_DV":JOINT_COLOR,"neither":FIG_STYLE["colors"]["other"]}

# Figure S7 combines the current fixed null splits and the six matched family restrictions.
null_results = pd.read_csv(R / f"sgcb_{RECORD_ID}_readcount_null_metrics.csv")
null_de = null_results.loc[null_results.target.eq("Mean expression")].sort_values("split")
null_dv = null_results.loc[null_results.target.ne("Mean expression")].sort_values("split")
shape = population_family_scenarios.loc[population_family_scenarios.panel.eq("Positive") &
    population_family_scenarios.scenario.eq("shape")]
shape_calls = shape.pivot(index="design", columns="family", values="FP").reindex(CONFIG["designs"])[
    ["GG", "Gamma", "Lognormal"]]
cost = population_family_summary.copy()
cost["position"] = cost.panel.map({"Positive":0,"Observed zeros":3})+cost.design.map(
    {"10:10":0,"10:50":1,"50:50":2})
cost["family_offset"] = cost.family.map({"GG":-.2,"Gamma":0,"Lognormal":.2})
cost["color"] = cost.family.map({"GG":DV_COLOR,"Gamma":NORMAL,"Lognormal":JOINT_COLOR})
fig, axes = plt.subplots(2,2,figsize=(7.25,5.9),layout="constrained")
axes[0,0].plot(null_de.split, null_de.raw_p005_full, "o", color=MEAN_COLOR, label="Mean")
axes[0,0].plot(null_dv.split, null_dv.raw_p005_full, "s", color=DV_COLOR, label="CV²")
axes[0,0].axhline(.05,color=RED,ls=":",lw=.8)
axes[0,0].set(title="$\\mathbf{a}$  Fixed GTEx liver nulls",xlabel="Partition",ylabel="Fraction with raw p < 0.05",
    xticks=np.arange(1,11),ylim=(0,max(.06,null_results.raw_p005_full.max()*1.15)))
axes[0,0].legend(frameon=False,fontsize=7)
axes[0,1].scatter(null_de.split-.08,null_de.calls,color=MEAN_COLOR,s=30,edgecolors="white",linewidths=.5,label="Mean")
axes[0,1].scatter(null_dv.split+.08,null_dv.calls,color=DV_COLOR,marker="s",s=30,edgecolors="white",linewidths=.5,label="CV²")
axes[0,1].set(title="$\\mathbf{b}$  Full-family BH discoveries",xlabel="Partition",ylabel="Genes called",
    xticks=np.arange(1,11),ylim=(-.5,max(1,null_results.calls.max()*1.15)))
axes[0,1].text(.97,.95,
    f"Partitions with ≥1 call\nMean: {int(null_de.calls.gt(0).sum())}/{len(null_de)} · CV²: {int(null_dv.calls.gt(0).sum())}/{len(null_dv)}",
    transform=axes[0,1].transAxes,ha="right",va="top",fontsize=7.5)
axes[1,0].scatter(np.repeat(np.arange(3),3)+np.tile([-.2,0,.2],3),shape_calls.to_numpy().ravel(),
    c=np.tile([DV_COLOR,NORMAL,JOINT_COLOR],3),s=28)
axes[1,0].set(title="$\\mathbf{c}$  Equal-mean, equal-CV² shape change",xlabel="Samples per group",
    ylabel="False calls / 600 genes",xticks=np.arange(3),xticklabels=CONFIG["designs"],
    ylim=(-.5,max(1,shape_calls.max().max()*1.15)))
axes[1,1].scatter(cost.position+cost.family_offset,cost.seconds,c=cost.color,s=28)
axes[1,1].set(title="$\\mathbf{d}$  DV component wall time",xlabel="Panel / sample design",
    ylabel="DV fitting wall time (s)",xticks=np.arange(6),
    xticklabels=["Positive\n10:10","Positive\n10:50","Positive\n50:50","Zeros\n10:10","Zeros\n10:50","Zeros\n50:50"])
axes[1,1].tick_params(axis="x",labelsize=7.2)
axes[0,0].grid(axis="y")
axes[0,1].grid(axis="y")
axes[1,0].grid(axis="y")
axes[1,1].grid(axis="y")
fig.legend(handles=[Line2D([],[],color=DV_COLOR,marker="o",ls="",label="GG"),
    Line2D([],[],color=NORMAL,marker="o",ls="",label="Gamma"),
    Line2D([],[],color=JOINT_COLOR,marker="o",ls="",label="Lognormal")],
    frameon=False,ncols=3,loc="outside lower center",fontsize=8)
fig.savefig(F / "Figure_S7.pdf")
fig.savefig(F / "Figure_S7.svg")
fig.savefig(F / "Figure_S7.png")
plt.close(fig)
print("Wrote Figures 1–2, 4–5 and S1–S7 as editable SVG/PDF and 300-dpi PNG.")
