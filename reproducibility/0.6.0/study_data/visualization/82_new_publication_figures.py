from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_rgba_array
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Patch

CONFIG = {"paths": {"root": Path(sys.argv[1])}, "threads": int(os.environ["SGCB_THREADS"]), "plot": {
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
    "savefig.dpi": 300, "figure.constrained_layout.h_pad": 5 / 72,
    "figure.constrained_layout.w_pad": 4 / 72},
    "figure_palettes": {'3': {'name': 'Sage and apricot',
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
    "distribution": {"bandwidth_adjust": .85, "width": .70, "point_size": 1.4},
    "heatmap": {"cell_gap": 1.8, "font_size": 7.2, "colorbar_height": .025,
                "panel_title_y":1.08,
                "truth_band": [1.012,0,.018,1], "truth_colors": ["#D5DDD6","#D7BF8A","#7C97AA"]},
    "zero_figure": {"figsize": [7.25,7.8], "height_ratios": [1.10,.85,2.1],
                    "method_labels": ["Gene-wise GG","SGCB","Gamma","Lognormal","Empirical delta","Jackknife"]},
    "detail": {"ink": "#192E38", "outline": "#4D6975", "interval_width": 1.5,
               "interval_caps": 2.7, "legend_size": 7.8,
               "feature_point_sizes": [2.5,4.5,6,6], "feature_point_alpha": [.34,.72,.80,.80]},
    "designs": ["10:10", "10:50", "50:50"],
    "methods": ["Gene-wise GG", "EB-GG-paired-population", "EB-Gamma-paired-population", "EB-Lognormal-paired-population", "Empirical-delta", "Empirical-jackknife"],
    "method_labels": ["GW-GG", "GG", "Gamma", "LN", "Delta", "JK"],
    "scenario_labels": ["LN null: pi .95", "LN null: pi .80", "LN null: pi .50", "LN null: pi .20",
        "Gamma null", "GG null", "Mixture null", "Near moment boundary", "Undefined second moment",
        "Zero frequency change", "Positive CV2 change", "Reinforcing changes", "Equal CV2: LN",
        "Equal CV2: Gamma", "Equal CV2: extreme pi", "Equal CV2: cross-family",
        "Mean only", "Mean + CV2", "LN null: no zeros", "CV2 change: no zeros"],
    "effect_models": ["SGCB", "unadjusted", "clinical", "clinical_and_composition"],
    "effect_labels": ["SGCB 95% interval", "Moment", "+ covariates", "+ composition"],
    "main_examples": ["IL1R2", "EPHX1"],
    "example_genes": ["TNFSF13B","IL1RN","IL1R2","MTCO1P40","HBA1","HBA2","HBB","SIAH2","ABCG8","EPHX1"]}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
plt.rcParams.update(CONFIG["plot"])
ROOT = CONFIG["paths"]["root"]
RECORD_ID = os.environ.get("SGCB_STUDY_RECORD_ID", "0.5.9")
R, F = ROOT / "results", ROOT / "manuscript" / "figures"
FIG_STYLE=CONFIG["figure_palettes"]["3"]
C=FIG_STYLE["colors"]

zero = pd.concat([pd.read_csv(R / "zero_calibration_n10_10_summary.csv"),
    pd.read_csv(R / "zero_calibration_n10_50_summary.csv"),
    pd.read_csv(R / "zero_calibration_n50_50_summary.csv")], ignore_index=True)
zero_scenarios = pd.concat([pd.read_csv(R / "zero_calibration_n10_10_scenarios.csv"),
    pd.read_csv(R / "zero_calibration_n10_50_scenarios.csv"),
    pd.read_csv(R / "zero_calibration_n50_50_scenarios.csv")], ignore_index=True)
zero["method"] = zero.method.replace({"SGCB":"Gene-wise GG"})
zero_scenarios["method"] = zero_scenarios.method.replace({"SGCB":"Gene-wise GG"})
zero.to_csv(R / "publication_zero_calibration.csv", index=False, encoding="utf-8")
zero_scenarios.to_csv(R / "publication_zero_calibration_scenarios.csv", index=False, encoding="utf-8")
population = pd.read_csv(R / "publication_population_family_summary.csv")
population = population.loc[population.panel.eq("Observed zeros") &
    population.family.eq("GG") & population.truth_domain.eq("finite_CV2")].copy()
population_scenarios = pd.read_csv(R / "publication_population_family_scenarios.csv")
population_scenarios = population_scenarios.loc[
    population_scenarios.panel.eq("Observed zeros") & population_scenarios.family.eq("GG")].copy()
restricted_population = pd.read_csv(R / "publication_population_family_summary.csv")
restricted_population = restricted_population.loc[restricted_population.panel.eq("Observed zeros") &
    restricted_population.truth_domain.eq("finite_CV2") & restricted_population.family.isin(["Gamma","Lognormal"])]
finite = pd.concat([zero[zero.truth_domain.eq("finite_CV2")], population,
    restricted_population], ignore_index=True).set_index(["design", "method"])
calls = finite.reindex(pd.MultiIndex.from_product([CONFIG["designs"], CONFIG["methods"]]))
gg = finite.xs("Gene-wise GG", level="method").reindex(CONFIG["designs"])
pop = finite.xs("EB-GG-paired-population", level="method").reindex(CONFIG["designs"])
sc = zero_scenarios[zero_scenarios.method.eq("Gene-wise GG")]
scenario_order = pd.read_csv(R / "zero_calibration_truth.csv").scenario.drop_duplicates().tolist()
sc10 = sc[sc.design.eq("10:10")].set_index("scenario").reindex(scenario_order)
sc50 = sc[sc.design.eq("50:50")].set_index("scenario").reindex(scenario_order)
pop10 = population_scenarios[population_scenarios.design.eq("10:10")].set_index("scenario").reindex(scenario_order)
pop50 = population_scenarios[population_scenarios.design.eq("50:50")].set_index("scenario").reindex(scenario_order)
scenario_labels = CONFIG["scenario_labels"]

# Figure 3 uses the same complete comparisons, with method rows and explicit truth classes.
fig=plt.figure(figsize=CONFIG["zero_figure"]["figsize"],layout="constrained")
top,middle,bottom=fig.subfigures(3,1,height_ratios=CONFIG["zero_figure"]["height_ratios"])
top.suptitle(r"$\mathbf{a}$",x=.015,ha="left",fontsize=10)
a0,a1,a2=top.subplots(1,3,sharex=True,sharey=True)
r0=finite.xs("10:10",level="design").reindex(CONFIG["methods"])
r1=finite.xs("10:50",level="design").reindex(CONFIG["methods"])
r2=finite.xs("50:50",level="design").reindex(CONFIG["methods"])
method_y=np.arange(len(CONFIG["methods"]))
bar_total=calls.TP+calls.FP+calls.calls_with_undefined_truth.fillna(0)
bar_limit=bar_total.max()*1.30
a0.barh(method_y,r0.TP,color=C["DE"],height=.66)
a0.barh(method_y,r0.FP,left=r0.TP,color=C["false"],height=.66)
a0.barh(method_y,r0.calls_with_undefined_truth.fillna(0),left=r0.TP+r0.FP,
    color=C["other"],edgecolor=C["control"],hatch="///",linewidth=.4,height=.66)
a1.barh(method_y,r1.TP,color=C["DE"],height=.66)
a1.barh(method_y,r1.FP,left=r1.TP,color=C["false"],height=.66)
a1.barh(method_y,r1.calls_with_undefined_truth.fillna(0),left=r1.TP+r1.FP,
    color=C["other"],edgecolor=C["control"],hatch="///",linewidth=.4,height=.66)
a2.barh(method_y,r2.TP,color=C["DE"],height=.66)
a2.barh(method_y,r2.FP,left=r2.TP,color=C["false"],height=.66)
a2.barh(method_y,r2.calls_with_undefined_truth.fillna(0),left=r2.TP+r2.FP,
    color=C["other"],edgecolor=C["control"],hatch="///",linewidth=.4,height=.66)
np.frompyfunc(a0.text,3,1)(r0.TP+r0.FP+r0.calls_with_undefined_truth.fillna(0)+bar_limit*.018,
    method_y,np.char.add(np.char.add(np.char.mod("%d",r0.TP.to_numpy()),"/"),np.char.mod("%d",r0.FP.to_numpy())))
np.frompyfunc(a1.text,3,1)(r1.TP+r1.FP+r1.calls_with_undefined_truth.fillna(0)+bar_limit*.018,
    method_y,np.char.add(np.char.add(np.char.mod("%d",r1.TP.to_numpy()),"/"),np.char.mod("%d",r1.FP.to_numpy())))
np.frompyfunc(a2.text,3,1)(r2.TP+r2.FP+r2.calls_with_undefined_truth.fillna(0)+bar_limit*.018,
    method_y,np.char.add(np.char.add(np.char.mod("%d",r2.TP.to_numpy()),"/"),np.char.mod("%d",r2.FP.to_numpy())))
plt.setp(list(a0.texts)+list(a1.texts)+list(a2.texts),fontsize=6.9,va="center",ha="left",color=CONFIG["detail"]["ink"])
plt.setp([a0,a1,a2],xlim=(0,bar_limit),ylim=(len(method_y)-.45,-.55),
    yticks=method_y,xlabel="Genes called")
a0.set(title="10 : 10",yticklabels=CONFIG["zero_figure"]["method_labels"])
a1.set(title="10 : 50")
a2.set(title="50 : 50")
a0.get_yticklabels()[1].set_fontweight("bold")
a0.tick_params(axis="y",labelsize=7.3,length=0)
a1.tick_params(axis="y",labelleft=False,length=0)
a2.tick_params(axis="y",labelleft=False,length=0)
a0.xaxis.set_major_locator(MaxNLocator(3,integer=True))
a1.xaxis.set_major_locator(MaxNLocator(3,integer=True))
a2.xaxis.set_major_locator(MaxNLocator(3,integer=True))
a0.grid(axis="x")
a1.grid(axis="x")
a2.grid(axis="x")
top.legend(handles=[Patch(facecolor=C["DE"],label="True positive"),
    Patch(facecolor=C["false"],label="False positive"),
    Patch(facecolor=C["other"],edgecolor=C["control"],hatch="///",label="Undefined CV² truth")],
    loc="outside lower center",ncols=3,frameon=False,fontsize=7.2)

middle.suptitle(r"$\mathbf{b}$",x=.015,ha="left",fontsize=10)
b0,b1,b2=middle.subplots(1,3,sharey=True)
b0.plot(np.arange(3),gg.coverage_available,"o-",color=C["DE"],ms=4)
b0.plot(np.arange(3),pop.coverage_available,"D-",color=C["DV"],ms=4)
b1.plot(np.arange(3),gg.coverage_all,"o-",color=C["DE"],ms=4)
b1.plot(np.arange(3),pop.coverage_all,"D-",color=C["DV"],ms=4)
b2.plot(np.arange(3),gg.availability,"o-",color=C["DE"],ms=4)
b2.plot(np.arange(3),pop.availability,"D-",color=C["DV"],ms=4)
b0.axhline(.95,color=C["control"],ls=":",lw=.75)
plt.setp([b0,b1,b2],xlim=(-.10,2.10),ylim=(.65,1.02),xticks=np.arange(3),
    xticklabels=CONFIG["designs"],yticks=[.7,.85,1],xlabel="Samples per group")
b0.set(ylabel="Covered /\navailable")
b1.set(ylabel="Covered /\nall finite")
b2.set(ylabel="Available /\nall finite")
b0.grid(True)
b1.grid(True)
b2.grid(True)
b0.legend(handles=[Line2D([],[],color=C["DE"],marker="o",ms=4,label="Gene-wise GG"),
    Line2D([],[],color=C["DV"],marker="D",ms=4,label="SGCB")],
    loc="lower left",ncols=1,frameon=False,fontsize=7.0,handlelength=1.3,labelspacing=.2,
    title="3,800 finite truths",title_fontsize=6.8)

heat_grid=bottom.add_gridspec(2,2,height_ratios=[1,CONFIG["heatmap"]["colorbar_height"]],
    width_ratios=[1.45,1])
c=bottom.add_subplot(heat_grid[0,0])
d=bottom.add_subplot(heat_grid[0,1])
availability_bar=bottom.add_subplot(heat_grid[1,0])
calls_bar=bottom.add_subplot(heat_grid[1,1])
availability=np.column_stack([sc10.available/sc10.features,sc50.available/sc50.features,
    pop10.available/pop10.features,pop50.available/pop50.features])
sns.heatmap(availability,ax=c,vmin=0,vmax=1,
    cmap=sns.blend_palette(FIG_STYLE["availability"],as_cmap=True),
    annot=np.char.mod("%.0f",availability*100),fmt="",annot_kws={"fontsize":CONFIG["heatmap"]["font_size"]},
    linewidths=CONFIG["heatmap"]["cell_gap"],linecolor="white",cbar_ax=availability_bar,
    cbar_kws={"orientation":"horizontal","ticks":[0,.5,1]},
    xticklabels=["10:10","50:50","10:10","50:50"],yticklabels=scenario_labels)
c.set(title=r"$\mathbf{c}$",xlabel="",ylabel="")
c.tick_params(axis="y",labelsize=7.0,length=0,labelrotation=0)
c.tick_params(axis="x",labelsize=7.3,length=0,labelrotation=0)
availability_header=c.secondary_xaxis("top")
availability_header.set_xticks([1,3],labels=["Gene-wise GG","SGCB"])
availability_header.tick_params(length=0,labelsize=7.3,pad=4)
availability_header.spines["top"].set_visible(False)
c.axvline(2,color="white",lw=2.2)
truth_summary=pd.read_csv(R/"zero_calibration_truth.csv").drop_duplicates("scenario").set_index("scenario").reindex(scenario_order)
truth_kind=np.select([~truth_summary.finite_second_moment.to_numpy(dtype=bool),
    truth_summary.dv_truth.astype(str).str.casefold().eq("true")],[2,1],default=0)
truth_band=c.inset_axes(CONFIG["heatmap"]["truth_band"])
truth_band.set_in_layout(False)
sns.heatmap(truth_kind[:,None],ax=truth_band,vmin=0,vmax=2,
    cmap=ListedColormap(CONFIG["heatmap"]["truth_colors"]),cbar=False,
    xticklabels=False,yticklabels=False,linewidths=.35,linecolor="white")
current_zero_scenarios=pd.read_csv(R/"publication_current_api_scenarios.csv")
current_zero_scenarios=current_zero_scenarios.loc[current_zero_scenarios.panel.eq("zero") & current_zero_scenarios.design.eq("50:50")].copy()
current_zero_scenarios["call_fraction"]=current_zero_scenarios.calls/current_zero_scenarios.features
zero_call_matrix=current_zero_scenarios.pivot(index="scenario",columns="endpoint",values="call_fraction").reindex(scenario_order)[["DE","DV"]]
call_labels=np.where(zero_call_matrix.to_numpy()==0,"0",np.char.mod("%.1f",100*zero_call_matrix.to_numpy()))
sns.heatmap(zero_call_matrix,ax=d,vmin=0,vmax=1,cmap=sns.blend_palette(FIG_STYLE["calls"],as_cmap=True),
    annot=call_labels,fmt="",annot_kws={"fontsize":CONFIG["heatmap"]["font_size"]},
    linewidths=CONFIG["heatmap"]["cell_gap"],linecolor="white",cbar_ax=calls_bar,
    cbar_kws={"orientation":"horizontal","ticks":[0,.5,1]},
    xticklabels=["Mean","CV²"],yticklabels=False)
d.set(title=r"$\mathbf{d}$",xlabel="50:50; 200 genes per scenario",ylabel="")
d.tick_params(length=0,labelrotation=0)
c.grid(False)
d.grid(False)
availability_bar.set_xticklabels(["0","50","100%"])
calls_bar.set_xticklabels(["0","50","100%"])
availability_bar.tick_params(length=2,pad=2,labelsize=7)
calls_bar.tick_params(length=2,pad=2,labelsize=7)
availability_bar.set_xlabel("Available comparisons (%)",fontsize=7.3,labelpad=3)
calls_bar.set_xlabel("Called genes (%)",fontsize=7.3,labelpad=3)
c.set_title(r"$\mathbf{c}$",y=CONFIG["heatmap"]["panel_title_y"])
d.set_title(r"$\mathbf{d}$",y=CONFIG["heatmap"]["panel_title_y"])
c.hlines([9,12,16,18],0,4,color="white",lw=1.8)
d.hlines([9,12,16,18],0,2,color="white",lw=1.8)
bottom.legend(handles=[Patch(facecolor=CONFIG["heatmap"]["truth_colors"][0],label="CV² null"),
    Patch(facecolor=CONFIG["heatmap"]["truth_colors"][1],label="CV² change"),
    Patch(facecolor=CONFIG["heatmap"]["truth_colors"][2],label="Undefined CV²")],
    title="Truth stripe",title_fontsize=7,loc="outside lower center",ncols=3,
    frameon=False,fontsize=7,labelspacing=.3)
fig.savefig(F / "Figure_3.svg")
fig.savefig(F / "Figure_3.pdf")
fig.savefig(F / "Figure_3.png",dpi=300)
plt.close(fig)

FIG_STYLE=CONFIG["figure_palettes"]["6"]
C=FIG_STYLE["colors"]
COHORT_COLORS={0:FIG_STYLE["cohorts"]["normal"],1:FIG_STYLE["cohorts"]["tumor"]}

# The fixed independent disease examples use the same public API and all profiles.
sle = pd.read_csv(R / f"sgcb_{RECORD_ID}_sle_classical_monocytes_2021_biological_results.csv")
liver = pd.read_csv(R / f"sgcb_{RECORD_ID}_masld_liver_2023_biological_results.csv")
sle = sle.loc[np.isfinite(sle.population_effect) & sle.population_available].copy()
liver = liver.loc[np.isfinite(liver.population_effect) & liver.population_available].copy()
sle["display_class"] = np.select([sle.population_DV & sle.DE,sle.population_DV,sle.DE],[3,2,1],default=0)
liver["display_class"] = np.select([liver.population_DV & liver.DE,liver.population_DV,liver.DE],[3,2,1],default=0)
sle = sle.sort_values("display_class")
liver = liver.sort_values("display_class")
palette = np.array([C["other"],C["DE"],C["DV"],C["both"]])
obs_sle = pd.read_csv(R / f"sgcb_{RECORD_ID}_sle_classical_monocytes_2021_example_observations.csv")
obs_liver = pd.read_csv(R / f"sgcb_{RECORD_ID}_masld_liver_2023_example_observations.csv")
fx_sle = pd.read_csv(R / f"sgcb_{RECORD_ID}_sle_classical_monocytes_2021_example_effects.csv")
fx_liver = pd.read_csv(R / f"sgcb_{RECORD_ID}_masld_liver_2023_example_effects.csv")
fig,ax = plt.subplots(3,2,figsize=(7.25,7.2),layout="constrained",
    gridspec_kw={"height_ratios":[1.1,1.05,.70]})
sle_colors=to_rgba_array(palette[sle.display_class])
liver_colors=to_rgba_array(palette[liver.display_class])
sle_colors[:,3]=np.take(CONFIG["detail"]["feature_point_alpha"],sle.display_class)
liver_colors[:,3]=np.take(CONFIG["detail"]["feature_point_alpha"],liver.display_class)
ax[0,0].scatter(sle.log2FoldChange,sle.population_effect,
    s=np.take(CONFIG["detail"]["feature_point_sizes"],sle.display_class),
    c=sle_colors,linewidths=0,rasterized=True)
ax[0,1].scatter(liver.log2FoldChange,liver.population_effect,
    s=np.take(CONFIG["detail"]["feature_point_sizes"],liver.display_class),
    c=liver_colors,linewidths=0,rasterized=True)
ax[0,0].set(title="$\\mathbf{a}$",
    xlabel="SGCB log₂ arithmetic-mean ratio",ylabel="SGCB log₂ CV² ratio",yscale="symlog")
ax[0,1].set(title="$\\mathbf{b}$",
    xlabel="SGCB log₂ arithmetic-mean ratio",ylabel="SGCB log₂ CV² ratio",yscale="symlog")
ax[0,0].axhline(0,color=C["control"],lw=.5,zorder=0)
ax[0,1].axhline(0,color=C["control"],lw=.5,zorder=0)
ax[0,0].text(.035,.96,f"SLE · classical monocytes\n24 controls / 64 patients\nShown: {len(sle):,} genes\nDV-only: {int(sle.display_class.eq(2).sum()):,}",
    transform=ax[0,0].transAxes,va="top",fontsize=7.2,
    bbox={"facecolor":"white","edgecolor":"none","alpha":.88,"pad":2})
ax[0,1].text(.035,.96,f"MASLD · liver biopsies\n69 controls / 299 patients\nShown: {len(liver):,} genes\nDV-only: {int(liver.display_class.eq(2).sum()):,}",
    transform=ax[0,1].transAxes,va="top",fontsize=7.2,
    bbox={"facecolor":"white","edgecolor":"none","alpha":.88,"pad":2})
# Gene examples are shown in genome-wide context (Dreamlet Fig. 4F).
hero_sle=sle.loc[sle.gene_name.eq("IL1R2")].iloc[0]
hero_liver=liver.loc[liver.gene_name.eq("EPHX1")].iloc[0]
ax[0,0].scatter(hero_sle.log2FoldChange,hero_sle.population_effect,s=52,
    facecolors="white",edgecolors=C["DV"],linewidths=1.3,zorder=6)
ax[0,1].scatter(hero_liver.log2FoldChange,hero_liver.population_effect,s=52,
    facecolors="white",edgecolors=C["DV"],linewidths=1.3,zorder=6)
ax[0,0].annotate("IL1R2",(hero_sle.log2FoldChange,hero_sle.population_effect),xytext=(.97,.96),
    textcoords="axes fraction",ha="right",va="top",color=C["DV"],fontsize=8.5,weight="bold",arrowprops={"arrowstyle":"-","color":C["DV"],"lw":.7})
ax[0,1].annotate("EPHX1",(hero_liver.log2FoldChange,hero_liver.population_effect),xytext=(18,-26),
    textcoords="offset points",color=C["DV"],fontsize=8.5,weight="bold",arrowprops={"arrowstyle":"-","color":C["DV"],"lw":.7})

fig.legend(handles=[Line2D([],[],marker="o",color="none",markerfacecolor=C["other"],markeredgecolor="white",label="Neither"),
    Line2D([],[],marker="o",color="none",markerfacecolor=C["DE"],markeredgecolor="white",label="DE only"),
    Line2D([],[],marker="o",color="none",markerfacecolor=C["DV"],markeredgecolor="white",label="DV only"),
    Line2D([],[],marker="o",color="none",markerfacecolor=C["both"],markeredgecolor="white",label="Both")],
    loc="outside upper center",ncol=4,frameon=False,fontsize=CONFIG["detail"]["legend_size"])
main_sle = CONFIG["main_examples"][0]
main_liver = CONFIG["main_examples"][1]
main_obs_sle=obs_sle.loc[obs_sle.gene_name.eq(main_sle)]
main_obs_liver=obs_liver.loc[obs_liver.gene_name.eq(main_liver)]
sns.violinplot(data=main_obs_sle,x="log2_abundance",y="group",hue="group",order=[0,1],hue_order=[0,1],
    orient="h",palette=COHORT_COLORS,legend=False,inner="box",inner_kws={"box_width":3,"whis_width":.8},cut=0,density_norm="width",saturation=1,
    width=CONFIG["distribution"]["width"],bw_adjust=CONFIG["distribution"]["bandwidth_adjust"],linewidth=.95,linecolor=CONFIG["detail"]["outline"],ax=ax[1,0])
sns.violinplot(data=main_obs_liver,x="log2_abundance",y="group",hue="group",order=[0,1],hue_order=[0,1],
    orient="h",palette=COHORT_COLORS,legend=False,inner="box",inner_kws={"box_width":3,"whis_width":.8},cut=0,density_norm="width",saturation=1,
    width=CONFIG["distribution"]["width"],bw_adjust=CONFIG["distribution"]["bandwidth_adjust"],linewidth=.95,linecolor=CONFIG["detail"]["outline"],ax=ax[1,1])
sns.swarmplot(data=main_obs_sle,x="log2_abundance",y="group",order=[0,1],orient="h",color="#526F75",size=CONFIG["distribution"]["point_size"],alpha=.55,ax=ax[1,0])
sns.swarmplot(data=main_obs_liver,x="log2_abundance",y="group",order=[0,1],orient="h",color="#526F75",size=CONFIG["distribution"]["point_size"],alpha=.50,ax=ax[1,1])
ax[1,0].set(title=r"$\mathbf{c}$",xlabel="Log₂(normalized count + 0.5)",ylabel="",ylim=(1.55,-.9),yticks=[0,1],yticklabels=["Control","SLE"])
ax[1,1].set(title=r"$\mathbf{d}$",xlabel="Log₂(FPKM + 0.5)",ylabel="",ylim=(1.55,-.9),yticks=[0,1],yticklabels=["Control","MASLD"])
ax[1,0].text(.97,.97,main_sle,transform=ax[1,0].transAxes,ha="right",va="top",
    fontsize=8,fontweight="bold",color=C["DV"])
ax[1,1].text(.97,.97,main_liver,transform=ax[1,1].transAxes,ha="right",va="top",
    fontsize=8,fontweight="bold",color=C["DV"])
ax[1,0].text(.03,.96,f"CV² q {hero_sle.population_BH:.3g} · DV-only",transform=ax[1,0].transAxes,va="top",fontsize=8,color=CONFIG["detail"]["ink"],fontweight="semibold")
ax[1,0].text(.03,.86,f"Mean q {hero_sle.padj:.3g}",transform=ax[1,0].transAxes,va="top",fontsize=8,color=CONFIG["detail"]["ink"])
ax[1,1].text(.03,.96,f"CV² q {hero_liver.population_BH:.3g} · DV-only",transform=ax[1,1].transAxes,va="top",fontsize=8,color=CONFIG["detail"]["ink"],fontweight="semibold")
ax[1,1].text(.03,.86,f"Mean q {hero_liver.padj:.3g}",transform=ax[1,1].transAxes,va="top",fontsize=8,color=CONFIG["detail"]["ink"])
e0 = fx_sle.loc[fx_sle.gene_name.eq(main_sle)].set_index("model").reindex(CONFIG["effect_models"][:3])
e1 = fx_liver.loc[fx_liver.gene_name.eq(main_liver)].set_index("model").reindex(CONFIG["effect_models"])
ax[2,0].errorbar(e0.effect,np.arange(3),xerr=np.vstack([e0.effect-e0.lower,e0.upper-e0.effect]),
    fmt="o",color=C["DV"],elinewidth=CONFIG["detail"]["interval_width"],
    capsize=CONFIG["detail"]["interval_caps"],ms=4.8,mec="white",mew=.5)
ax[2,1].errorbar(e1.effect,np.arange(4),xerr=np.vstack([e1.effect-e1.lower,e1.upper-e1.effect]),
    fmt="o",color=C["DV"],elinewidth=CONFIG["detail"]["interval_width"],
    capsize=CONFIG["detail"]["interval_caps"],ms=4.8,mec="white",mew=.5)
ax[2,0].scatter(e0.effect.iloc[0],0,marker="D",s=34,color=C["DV"],edgecolors="white",linewidths=.5,zorder=4)
ax[2,1].scatter(e1.effect.iloc[0],0,marker="D",s=34,color=C["DV"],edgecolors="white",linewidths=.5,zorder=4)
ax[2,0].axhspan(-.35,.35,color=C["DV"],alpha=.055,zorder=0)
ax[2,1].axhspan(-.35,.35,color=C["DV"],alpha=.055,zorder=0)
ax[2,0].set(title="$\\mathbf{e}$",yticks=np.arange(3),
    yticklabels=CONFIG["effect_labels"][:3],xlabel=main_sle+" log₂ CV² ratio",ylim=(2.7,-.6))
ax[2,1].set(title="$\\mathbf{f}$",yticks=np.arange(4),
    yticklabels=CONFIG["effect_labels"],xlabel=main_liver+" log₂ CV² ratio",ylim=(3.7,-.6))
ax[2,0].axvline(0,color=C["control"],lw=.7,ls=":")
ax[2,1].axvline(0,color=C["control"],lw=.7,ls=":")
ax[0,0].grid(True)
ax[0,1].grid(True)
ax[1,0].grid(axis="x")
ax[1,1].grid(axis="x")
ax[2,0].grid(axis="x")
ax[2,1].grid(axis="x")
fig.savefig(F / "Figure_6.svg")
fig.savefig(F / "Figure_6.pdf")
fig.savefig(F / "Figure_6.png",dpi=300)
plt.close(fig)

# All ten prespecified examples remain visible, including genes without a DV call.
fig,ax = plt.subplots(5,2,figsize=(7.25,10.1),layout="constrained",sharey=True)
example = CONFIG["example_genes"][0]
ax[0,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[0,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[0,0].set(title="$\\mathbf{a}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][1]
ax[0,1].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[0,1].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[0,1].set(title="$\\mathbf{b}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][2]
ax[1,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[1,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[1,0].set(title="$\\mathbf{c}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][3]
ax[1,1].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[1,1].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[1,1].set(title="$\\mathbf{d}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][4]
ax[2,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[2,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[2,0].set(title="$\\mathbf{e}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][5]
ax[2,1].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[2,1].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[2,1].set(title="$\\mathbf{f}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][6]
ax[3,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[3,0].ecdf(obs_sle.loc[obs_sle.gene_name.eq(example) & obs_sle.group.eq(1),"log2_abundance"],
    color=C["case"],label="SLE")
ax[3,0].set(title="$\\mathbf{g}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][7]
ax[3,1].ecdf(obs_liver.loc[obs_liver.gene_name.eq(example) & obs_liver.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[3,1].ecdf(obs_liver.loc[obs_liver.gene_name.eq(example) & obs_liver.group.eq(1),"log2_abundance"],
    color=C["case"],linestyle="--",label="MASLD")
ax[3,1].set(title="$\\mathbf{h}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][8]
ax[4,0].ecdf(obs_liver.loc[obs_liver.gene_name.eq(example) & obs_liver.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[4,0].ecdf(obs_liver.loc[obs_liver.gene_name.eq(example) & obs_liver.group.eq(1),"log2_abundance"],
    color=C["case"],linestyle="--",label="MASLD")
ax[4,0].set(title="$\\mathbf{i}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
example = CONFIG["example_genes"][9]
ax[4,1].ecdf(obs_liver.loc[obs_liver.gene_name.eq(example) & obs_liver.group.eq(0),"log2_abundance"],
    color=C["control"],label="Control")
ax[4,1].ecdf(obs_liver.loc[obs_liver.gene_name.eq(example) & obs_liver.group.eq(1),"log2_abundance"],
    color=C["case"],linestyle="--",label="MASLD")
ax[4,1].set(title="$\\mathbf{j}$  "+example,xlabel="Log₂(abundance + 0.5)",
    ylabel="Cumulative fraction",ylim=(0,1.02))
plt.setp(ax.ravel(),xlabel="",ylabel="",yticks=[0,.5,1])
fig.supxlabel("Log₂(abundance + 0.5)",fontsize=9.2)
fig.supylabel("Cumulative fraction",fontsize=9.2)
fig.legend(handles=[Line2D([],[],color=C["control"],label="Control"),
    Line2D([],[],color=C["case"],label="SLE"),
    Line2D([],[],color=C["case"],ls="--",label="MASLD")],
    frameon=False,ncols=3,loc="outside upper center",fontsize=CONFIG["detail"]["legend_size"])
fig.savefig(F / "Figure_S8.svg")
fig.savefig(F / "Figure_S8.pdf")
fig.savefig(F / "Figure_S8.png",dpi=300)
plt.close(fig)
print("Exported observed-zero calibration, current disease applications and all prespecified examples")
