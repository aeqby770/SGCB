from pathlib import Path
import os
import numpy as np
import pandas as pd
root=Path(os.environ["SGCB_R2_ROOT"]);r=root/"results"
def load(name):return pd.read_csv(r/name,sep="\t",keep_default_na=False,na_values=[""])
def pareto(scope,objectives,signs):
    arr=objectives.to_numpy()*np.asarray(signs)
    assert np.isfinite(arr).all()
    ge=(arr[:,None,:]>=arr[None,:,:]).all(axis=2)
    gt=(arr[:,None,:]>arr[None,:,:]).any(axis=2)
    dom = ge & gt
    return pd.DataFrame({"scope":scope,"method":objectives.index,"n_objectives":arr.shape[1],
        "nondominated_point_estimate":~dom.any(axis=0),
        "dominators":[", ".join(objectives.index[dom[:,j]]) for j in range(len(objectives))]})
def write_matrix(name,df):
    df.rename_axis("method").reset_index().to_csv(r/(name+".tsv"),sep="\t",index=False)
h=load("mean_candidate_holdout.tsv"); parts=[]
for n in ["10_10","10_50","50_50"]:
    for panel,what in [("mixed","power"),("mixed","FDP"),("null","any_call")]:
        a=h[(h["size"]==n)&(h.panel==panel)&(h.variable==what)].set_index("method").value
        parts.append(a.rename(n+"_"+panel+"_"+what))
mean=pd.concat(parts,axis=1)
write_matrix("pareto_mean_objectives",mean)
ans=[pareto("arithmetic mean, held-out independent matrices",mean,[1,-1,-1]*3)]
s=load("independent_operating_characteristics.tsv")
s=s[(s.state=="all")&(s.endpoint=="dv")&s.method.isin(["SGCB","Moments"])]
parts=[]
for n in ["10_10","10_50","50_50"]:
    for panel,what in [("mixed","power"),("mixed","FDP"),("mixed","interval_score"),("null","any_call")]:
        a=s[(s["size"]==n)&(s.panel==panel)&(s.metric==what)].set_index("method").value
        parts.append(a.rename(n+"_"+panel+"_"+what))
variability=pd.concat(parts,axis=1)
write_matrix("pareto_variability_objectives",variability)
ans.append(pareto("observed CV2, independent matrices",variability,[1,-1,-1,-1]*3))
o=load("outer_refit_summary.tsv")
cols={}
for name in ["conditional","outer"]:
    result={}
    for _,v in o.iterrows():
        k=v["size"]+"_"+v.state
        result[k+"_coverage_error"]=abs(v[name+"_coverage"]-.95)
        result[k+"_width"]=v[name+"_width"]
        result[k+"_interval_score"]=v[name+"_score"]
    cols[name]=result
intervals=pd.DataFrame.from_dict(cols,orient="index");write_matrix("pareto_interval_objectives",intervals)
ans.append(pareto("conditional intervals on the same probes",intervals,[-1]*len(intervals.columns)))
depth=pd.concat([load("exposure_summary.tsv").query("scope == 'all matrices'"),
    load("known_exposure_ql_summary.tsv")],ignore_index=True)
f=depth[depth.condition=="technical_only"].set_index("arm")
m=depth[depth.condition=="technical_and_biological"].set_index("arm")
mat=pd.concat([f.any_probability.rename("null_any"),f.false_calls.rename("null_false_call_count"),m.mean_FDP.rename("mixed_FDP"),m.false_calls.rename("mixed_false_calls"),m.mean_power.rename("mixed_power")],axis=1)
write_matrix("pareto_known_depth_objectives",mat)
ans.append(pareto("arithmetic means in exposure-imbalanced counts",mat,[-1,-1,-1,-1,1]))
boundary=load("boundary_operating_characteristics.tsv")
old_new={}
for nm in ["released","ratio_limit"]:
    z=boundary[boundary.candidate==nm].set_index(["stage","arm"])
    old_new[nm]={
        "fresh_null_any":z.loc[("fresh_depth","depth_null"),"family_probability"],
        "fresh_null_false_calls":z.loc[("fresh_depth","depth_null"),"false_calls"],
        "fresh_mixed_FDP":z.loc[("fresh_depth","depth_mixed"),"mean_FDP"],
        "fresh_mixed_power":z.loc[("fresh_depth","depth_mixed"),"power"],
        "all_twofold_power":z.loc[("positive","all_double"),"power"],
        "all_fourfold_power":z.loc[("positive","all_quadruple"),"power"],
        "balanced_fourfold_power":z.loc[("positive","balanced_fourfold"),"power"]}
matrix=pd.DataFrame.from_dict(old_new,orient="index")
write_matrix("pareto_version_objectives",matrix)
ans.append(pareto("R1 and 0.6.0 under technical and biological boundary stress",matrix,[-1,-1,-1,1,1,1,1]))
pd.concat(ans,ignore_index=True).to_csv(r/"empirical_pareto_sets.tsv",sep="\t",index=False)
print(pd.concat(ans,ignore_index=True).to_string(index=False))
