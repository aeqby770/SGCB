from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score

CONFIG = {"paths": {"root": Path(__file__).resolve().parents[1],
                   "mean_prefix": os.environ.get("SGCB_MEAN_PREFIX", "sgcb_0.5.7_seqc")},
          "threads": int(os.environ["SGCB_THREADS"]), "plot": {}, "thresholds": [0.5, 1.0, 1.5]}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
plt.rcParams.update(CONFIG["plot"])
ROOT = CONFIG["paths"]["root"]
OUTPUT = ROOT / "results" / (CONFIG["paths"]["mean_prefix"] + "_qpcr")
ref = pd.read_csv(ROOT / "results/taqman_reference_all_assays.csv")
old = pd.read_csv(ROOT / "data/seqc_archived_comparator_scores.csv", dtype={"gene": str})
# The archived glmGamPoi wrapper kept its sequential result names.
# Its baseMean vector and row order identify the shared count-matrix rows.
glm_rows = old.method.eq("glmGamPoi")
common_rows = old.loc[old.method.eq("edgeR"), "gene"].to_numpy()
assert np.array_equal(old.loc[glm_rows, "gene"].astype(int).to_numpy(),
                      np.arange(1, len(common_rows) + 1))
assert np.array_equal(old.loc[glm_rows, "baseMean"].to_numpy(),
                      old.loc[old.method.eq("edgeR"), "baseMean"].to_numpy())
old.loc[glm_rows, "gene"] = common_rows
new = pd.read_csv(ROOT / "results" / (CONFIG["paths"]["mean_prefix"] + "_genes.csv"))
sgcb = new[["gene_id", "pvalue", "padj", "log2FoldChange"]].rename(
    columns={"gene_id": "gene", "log2FoldChange": "log2FC"}).assign(method="SGCB")
scores = pd.concat([old.loc[~old.method.isin(["samr", "SAMseq"]),
                           ["gene", "pvalue", "padj", "log2FC", "method"]], sgcb])
wide_p = scores.pivot(index="gene", columns="method", values="pvalue").dropna()
wide_q = scores.pivot(index="gene", columns="method", values="padj").loc[wide_p.index]
wide_fc = scores.pivot(index="gene", columns="method", values="log2FC").loc[wide_p.index]
eligible = ref.loc[ref.all_eight_detected & ~ref.duplicate_symbol].set_index("gene")
genes = wide_p.index.intersection(eligible.index, sort=True)
reference = eligible.loc[genes]
methods = wide_p.columns.to_numpy()
score = -np.log10(np.clip(wide_p.loc[genes].to_numpy(), 1e-300, 1))
prediction = wide_q.loc[genes].to_numpy() < 0.05
reference_fc = np.column_stack([reference.log2FC_qpcr_B_over_A.to_numpy(),
    np.log2(reference[["B1", "B2", "B3", "B4"]].mean(axis=1).to_numpy() /
            reference[["A1", "A2", "A3", "A4"]].mean(axis=1).to_numpy())])
reference_names = np.array(["geometric_mean", "arithmetic_mean"])
reference_types = np.repeat(reference_names, len(CONFIG["thresholds"]))
thresholds = np.tile(np.asarray(CONFIG["thresholds"]), len(reference_names))
# Preserve the original geometric reference and add the matched arithmetic
# sensitivity on the identical genes, methods and three thresholds.
truth = np.repeat(np.abs(reference_fc), len(CONFIG["thresholds"]), axis=1) >= thresholds
n = len(genes)
n_pos = truth.sum(axis=0)
n_neg = (~truth).sum(axis=0)
# DeLong placements share exactly the same qPCR genes across methods.
comparison = (score[:, None, :] > score[None, :, :]).astype(float)
comparison += 0.5 * (score[:, None, :] == score[None, :, :])
pos_placement = np.einsum("ijm,jt->imt", comparison, (~truth).astype(float)) / n_neg
neg_placement = np.einsum("ijm,it->jmt", comparison, truth.astype(float)) / n_pos
auc = np.einsum("imt,it->mt", pos_placement, truth.astype(float)) / n_pos
var_pos = np.einsum("imt,it->mt", (pos_placement - auc[None, :, :]) ** 2,
                    truth.astype(float)) / (n_pos - 1)
var_neg = np.einsum("imt,it->mt", (neg_placement - auc[None, :, :]) ** 2,
                    (~truth).astype(float)) / (n_neg - 1)
auc_se = np.sqrt(var_pos / n_pos + var_neg / n_neg)
tp = (prediction[:, :, None] & truth[:, None, :]).sum(axis=0)
fp = (prediction[:, :, None] & ~truth[:, None, :]).sum(axis=0)
fn = (~prediction[:, :, None] & truth[:, None, :]).sum(axis=0)
tn = (~prediction[:, :, None] & ~truth[:, None, :]).sum(axis=0)
f1 = 2 * tp / (2 * tp + fp + fn)
prob = np.stack([tp, fp, fn, tn], axis=-1) / n
a, b, c = prob[:, :, 0], prob[:, :, 1], prob[:, :, 2]
denom = 2 * a + b + c
gradient = np.stack([2 * (b + c) / denom ** 2, -2 * a / denom ** 2,
                     -2 * a / denom ** 2, np.zeros_like(a)], axis=-1)
f1_se = np.sqrt(((gradient ** 2 * prob).sum(axis=-1) -
                 (gradient * prob).sum(axis=-1) ** 2) / n)
table = pd.DataFrame({"method": np.repeat(methods, len(thresholds)),
    "reference_estimand": np.tile(reference_types, len(methods)),
    "qpcr_abs_log2FC_threshold": np.tile(thresholds, len(methods)), "n_common_genes": n,
    "n_positive": np.tile(n_pos, len(methods)), "n_negative": np.tile(n_neg, len(methods)),
    "AUC": auc.ravel(), "AUC_low95": np.maximum(0, auc - 1.96 * auc_se).ravel(),
    "AUC_high95": np.minimum(1, auc + 1.96 * auc_se).ravel(),
    "TP": tp.ravel(), "FP": fp.ravel(), "FN": fn.ravel(), "TN": tn.ravel(),
    "F1": f1.ravel(), "F1_low95": np.maximum(0, f1 - 1.96 * f1_se).ravel(),
    "F1_high95": np.minimum(1, f1 + 1.96 * f1_se).ravel(),
    "precision": (tp / (tp + fp)).ravel(), "recall": (tp / (tp + fn)).ravel()})
table["average_precision"] = average_precision_score(
    np.broadcast_to(truth[:, None, :], (n, len(methods), len(thresholds))).reshape(n, -1),
    np.broadcast_to(score[:, :, None], (n, len(methods), len(thresholds))).reshape(n, -1),
    average=None)
table.to_csv(Path(str(OUTPUT) + "_metrics.csv"), index=False, encoding="utf-8")
merged = reference.reset_index().assign(
    log2FC_qpcr_arithmetic_B_over_A=reference_fc[:, 1]).merge(scores, on="gene", how="left")
merged.to_csv(Path(str(OUTPUT) + "_gene_scores.csv"), index=False, encoding="utf-8")
ref.assign(in_common_evaluation=ref.gene.isin(genes)).to_csv(
    Path(str(OUTPUT) + "_mapping.csv"), index=False, encoding="utf-8")
rho = spearmanr(reference_fc, wide_fc.loc[genes, "SGCB"]).statistic[:2, 2]
pd.DataFrame({"method": "SGCB", "contrast": "B_over_A",
    "reference_estimand": reference_names, "spearman_logFC": rho}).to_csv(
    Path(str(OUTPUT) + "_foldchange_correlations.csv"), index=False, encoding="utf-8")
print(table.loc[table.qpcr_abs_log2FC_threshold.eq(1)].to_string(index=False))
print({"SGCB_logFC_Spearman": rho})
print({"all_assays": len(ref), "all_eight_detected": int(ref.all_eight_detected.sum()),
       "unique_detected_assays": len(eligible), "common_evaluation": n})
