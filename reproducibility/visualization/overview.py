"""Render the SGCB mechanism figure for 30_publication_figures.py.

The caller supplies ROOT, FIG_STYLE and THREADS. The input/model/analysis
layout follows Dreamlet and DEVIL Fig. 1. All plotted glyphs are schematic.
"""
import os
from pathlib import Path as FilePath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, Rectangle, PathPatch, FancyArrowPatch, Polygon, Ellipse
from matplotlib.path import Path
from matplotlib.colors import ListedColormap, to_rgb
from scipy.special import gammaln, gammaincinv
import seaborn as sns

CONFIG={"paths":{"root":ROOT,"fonts":ROOT/"visualization"/"_fonts",
        "system_fonts":FilePath(os.environ["WINDIR"])/"Fonts"},
    "threads":THREADS,"plot":{"font.family":["Segoe Print","Inter","Arial"],
        "font.size":8.3,"text.color":"#20343E","axes.labelcolor":"#20343E",
        "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
        "mathtext.fontset":"stix","figure.facecolor":"#F8F4E9","savefig.facecolor":"#F8F4E9",
        "axes.facecolor":"#F8F4E9",
        "savefig.dpi":450},
    "canvas":{"width":7.25,"height":4.7},
    "layout":{"centers":[.67,1.98,3.50,5.14,6.63],"heading_y":4.08,
        "tube_x":[.24,.38,.52,.75,.89,1.03],"tube_y":3.20,
        "matrix":[.22,2.03,.94,.94],"law":[1.46,2.27,1.04,1.23],
        "parameter_icons":[1.46,1.28,1.04,.54],
        "population":[2.79,1.00,1.42,2.78],"state_grid":[2.91,2.03,1.18,1.12],
        "prior":[2.96,1.51,1.08,.27],
        "mean_icon":[4.57,2.77,1.18,.74],"cv_icon":[4.57,1.31,1.18,.74],
        "document":[6.12,1.12,.98,2.52],
        "links":[[[1.18,2.50],[1.43,2.50]],[[2.53,2.50],[2.76,2.50]],
                 [[4.24,2.50],[4.54,3.14]],[[4.24,2.50],[4.54,1.70]],
                 [[5.79,3.14],[6.08,3.14]],[[5.79,1.70],[6.08,1.70]]]},
    "style":{"ink":"#20343E","outline":"#52717E","muted":"#657D87",
        "light":"#E7EFF1","heading_size":9.0,"small_size":7.3,
        "arrow_width":1.05,"arrow_head":7.6,"tube_width":.105,"tube_height":.36,
        "tube_liquid":[.42,.58,.49,.57,.44,.53],
        "model_color":"#B58A43","location_color":"#627D94",
        "scale_color":"#6FABA1","shape_color":"#C79D58",
        "header_fill":"#E9F0E6","header_edge":"#D6E2D4",
        "depth":[.040,-.035],"board_front":"#F2F7ED",
        "board_side":"#BFCFBF","board_edge":"#A8BEAD",
        "paper_front":"#FFFEF7","paper_back":"#DCE7DE",
        "shadow_color":"#56706B"},
    "schematic":{"model_parameters":[1.1,.5,.7],
        "state_sigma":[.34,.5,.68],"state_Q":[.55,.8,1.1],
        "state_epsilon":[0,.05,.12],"state_relative_shift":[0,1.5,2.5],
        "state_nonzero":[.64,.84,.98],"band_probabilities":[.10,.25,.75,.90],
        "cloud":{"lower":.08,"upper":.92,"points":17},
        "weights":[[.03,.07,.03],[.06,.40,.11],[.02,.10,.18]],
        "contrast_shape":[9,9,9,2],"contrast_mean":[1.05,1.70,1.45,1.45]}}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
fm.fontManager.addfont(str(CONFIG["paths"]["fonts"] / "Inter-400.ttf"))
fm.fontManager.addfont(str(CONFIG["paths"]["fonts"] / "Inter-600.ttf"))
fm.fontManager.addfont(str(CONFIG["paths"]["fonts"] / "Montserrat-700.ttf"))
fm.fontManager.addfont(str(CONFIG["paths"]["system_fonts"] / "segoepr.ttf"))
fm.fontManager.addfont(str(CONFIG["paths"]["system_fonts"] / "segoeprb.ttf"))
previous_plot_style=plt.rcParams.copy()
plt.rcParams.update(CONFIG["plot"])
W,H=CONFIG["canvas"]["width"],CONFIG["canvas"]["height"]
INK=CONFIG["style"]["ink"]
OUTLINE=CONFIG["style"]["outline"]
MEAN,DV,CREAM=FIG_STYLE["colors"]["DE"],FIG_STYLE["colors"]["DV"],FIG_STYLE["colors"]["both"]
F=CONFIG["paths"]["root"] / "manuscript" / "figures"
fig=plt.figure(figsize=(W,H),dpi=450)
art=fig.add_axes([0,0,1,1])
art.set(xlim=(0,W),ylim=(0,H),aspect="equal")
art.axis("off")
art.add_patch(Rectangle((0,H-.43),W,.43,facecolor=CONFIG["style"]["header_fill"],
    edgecolor="none",zorder=-4))
art.hlines(H-.43,.14,W-.14,color=CONFIG["style"]["header_edge"],lw=.55,zorder=-3)

# Physical coordinates define the editable vector objects.
def fr(x,y,w,h):
    return [x/W,y/H,w/W,h/H]

def tube(x,color,liquid):
    y=CONFIG["layout"]["tube_y"]
    width,height=CONFIG["style"]["tube_width"],CONFIG["style"]["tube_height"]
    vertices=[(x,y+height),(x,y+.12),(x,y+.045),(x+width*.5,y),
        (x+width*.5,y),(x+width*.5,y),(x+width,y+.045),(x+width,y+.12),
        (x+width,y+height),(x,y+height)]
    codes=[Path.MOVETO,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,
        Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CLOSEPOLY]
    shape=Path(vertices,codes)
    shell=PathPatch(shape,facecolor="white",edgecolor=color,lw=.72,zorder=3)
    shell.set_path_effects([pe.SimplePatchShadow(offset=(1.2,-1.0),alpha=.12,
        shadow_rgbFace=CONFIG["style"]["shadow_color"]),pe.Normal()])
    art.add_patch(shell)
    liquid_patch=Rectangle((x,y),width,height*liquid,facecolor=color,
        edgecolor="none",alpha=.24,zorder=3.1)
    liquid_patch.set_clip_path(shell)
    art.add_patch(liquid_patch)
    liquid_core=Rectangle((x+width*.12,y),width*.76,height*liquid*.78,
        facecolor=color,edgecolor="none",alpha=.18,zorder=3.12)
    liquid_core.set_clip_path(shell)
    art.add_patch(liquid_core)
    art.plot([x+.012,x+width-.012],[y+height*liquid,y+height*liquid],
        color=color,lw=.55,alpha=.75,zorder=3.15)
    art.add_patch(PathPatch(shape,facecolor="none",edgecolor=color,lw=.72,zorder=3.2))
    meniscus=Ellipse((x+width/2,y+height*liquid),width*.86,.020,
        facecolor=color,edgecolor="white",lw=.28,alpha=.40,zorder=3.16)
    meniscus.set_clip_path(shell)
    art.add_patch(meniscus)
    art.add_patch(FancyBboxPatch((x-.016,y+height-.025),width+.032,.052,
        boxstyle="round,pad=0,rounding_size=.012",facecolor=color,edgecolor="white",lw=.3,zorder=4))
    cap_top=tuple(np.asarray(to_rgb(color))*.65+.35)
    art.add_patch(Ellipse((x+width/2,y+height+.026),width+.030,.021,
        facecolor=cap_top,edgecolor=color,lw=.35,zorder=4.15))
    art.plot([x+width*.27,x+width*.27],[y+.15,y+height-.07],
        color="white",alpha=.85,lw=.55,zorder=4)
    art.vlines(x+np.linspace(-.006,width+.006,5),y+height-.012,y+height+.018,
        color="white",alpha=.65,lw=.27,zorder=4.1)
    art.hlines(y+height*np.array([.36,.53,.70,.84]),x+width*.73,x+width*.94,
        color=color,alpha=.65,lw=.30,zorder=4.1)

def connector(start,end,color):
    dx,dy=end[0]-start[0],end[1]-start[1]
    path=Path([start,(start[0]+dx/2,start[1]+dy/6),
        (end[0]-dx/2,end[1]-dy/6),end],
        [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4])
    art.add_patch(FancyArrowPatch(path=path,arrowstyle="-|>",
        mutation_scale=CONFIG["style"]["arrow_head"],lw=CONFIG["style"]["arrow_width"],
        color=color,capstyle="round",joinstyle="round",zorder=6))

art.text(.16,4.48,"SGCB",fontfamily="Segoe Print",fontsize=17.5,fontweight="bold",va="center")
art.text(1.35,4.48,"Zero-aware mean and relative variability",fontsize=10.4,va="center")
headings=["Observed\nabundance","Zero-aware\ndistribution","Across-gene\nsharing",
    "Mean and\nvariability","Gene-level\nreadout"]
heading_text=np.frompyfunc(art.text,3,1)(np.array(CONFIG["layout"]["centers"]),
    np.repeat(CONFIG["layout"]["heading_y"],5),headings)
plt.setp(heading_text.tolist(),ha="center",va="center",fontsize=CONFIG["style"]["heading_size"],
    fontweight="semibold",linespacing=1.25)

# Input samples and matrix: the statistical unit and data object remain explicit.
tube_colors=[MEAN,MEAN,MEAN,DV,DV,DV]
art.add_patch(FancyBboxPatch((.205,3.105),.44,.067,
    boxstyle="round,pad=0,rounding_size=.012",facecolor="#B4D7D8",
    edgecolor="#78A8AD",lw=.5,zorder=2))
art.add_patch(Polygon([(.205,3.172),(.645,3.172),(.665,3.190),(.225,3.190)],
    facecolor="#E4F1EA",edgecolor="#89B2B0",lw=.45,zorder=2.1))
art.add_patch(FancyBboxPatch((.715,3.105),.44,.067,
    boxstyle="round,pad=0,rounding_size=.012",facecolor="#E7C1AF",
    edgecolor="#C89A84",lw=.5,zorder=2))
art.add_patch(Polygon([(.715,3.172),(1.155,3.172),(1.175,3.190),(.735,3.190)],
    facecolor="#F5E4D4",edgecolor="#D3AF94",lw=.45,zorder=2.1))
tuple(map(tube,CONFIG["layout"]["tube_x"],tube_colors,CONFIG["style"]["tube_liquid"]))
art.text(.42,3.66,"Control",ha="center",fontsize=7.5,color=MEAN,fontweight="semibold")
art.text(.92,3.66,"Case",ha="center",fontsize=7.5,color=DV,fontweight="semibold")
mx,my,mw,mh=CONFIG["layout"]["matrix"]
depth_x,depth_y=CONFIG["style"]["depth"]
art.add_patch(Polygon([(mx,my),(mx+mw,my),(mx+mw+depth_x,my+depth_y),(mx+depth_x,my+depth_y)],
    facecolor="#C2D7D0",edgecolor="#94B3AC",lw=.45,zorder=1.2))
art.add_patch(Polygon([(mx+mw,my),(mx+mw,my+mh),(mx+mw+depth_x,my+mh+depth_y),
    (mx+mw+depth_x,my+depth_y)],facecolor="#A7C2BD",edgecolor="#94B3AC",lw=.45,zorder=1.2))
input_ax=fig.add_axes(fr(mx,my,mw,mh))
gx,gy=np.meshgrid(np.arange(8),np.arange(8))
matrix=np.where((gx+2*gy)%7==0,0,np.where(gx<4,1+(gx+2*gy)%4,5+(gx+2*gy)%4))
matrix_colors=["#FBF7E9"]+sns.light_palette(MEAN,n_colors=5).as_hex()[1:]+sns.light_palette(DV,n_colors=5).as_hex()[1:]
sns.heatmap(matrix,ax=input_ax,
    cmap=ListedColormap(matrix_colors),vmin=0,vmax=8,
    linewidths=.65,linecolor="white",cbar=False,xticklabels=False,yticklabels=False)
input_ax.axvline(4,color="white",lw=1.25)
art.text(mx+mw/2,my+depth_y-.12,"Genes × samples",ha="center",fontsize=7.6,fontfamily="Inter")
art.plot([.37,.37,.95,.95],[3.10,3.05,3.05,3.10],color=OUTLINE,lw=.65)
connector((.66,3.05),(.66,2.99),OUTLINE)
art.text(.68,1.40,"TMM / input scale",ha="center",va="center",fontsize=7.4,fontfamily="Inter")
art.text(.68,.94,r"$X_{gj}=Y_{gj}/s_j$",ha="center",fontsize=10.2)

# One positive density explains the model; other glyphs use different objects.
law_ax=fig.add_axes(fr(*CONFIG["layout"]["law"]))
law_ax.set(xlim=(0,1),ylim=(0,1))
law_ax.axis("off")
support=np.linspace(.02,13,300)[:,None]
mu,sigma,Q=CONFIG["schematic"]["model_parameters"]
k=Q**-2
z=(np.log(support)-mu)/sigma
density=np.exp(np.log(Q)+k*np.log(k)-gammaln(k)-k-(np.expm1(Q*z)-Q*z)/Q**2)/(support*sigma)
curve_x=.23+.054*support[:,0]
curve_y=.16+.62*density/density.max()
law_ax.plot([.025,.98],[.16,.16],color=OUTLINE,lw=.65)
law_ax.plot([.10,.10],[.16,.68],color=OUTLINE,lw=1.7)
law_ax.scatter(.10,.68,color=OUTLINE,s=17,zorder=3)
law_ax.text(.10,.83,r"$1-\pi$",ha="center",fontsize=10)
law_ax.fill_between(curve_x,.16,curve_y[:,0],color=CONFIG["style"]["model_color"],alpha=.20)
law_ax.plot(curve_x,curve_y[:,0],color=CONFIG["style"]["model_color"],lw=1.45)
law_ax.text(.62,.91,r"$\pi f_+(x)$",ha="center",fontsize=10.3)
law_ax.text(.10,.02,"0",ha="center",fontsize=7.3,fontfamily="Inter")
law_ax.text(.66,.02,"Positive",ha="center",fontsize=7.3)
profile_ax=fig.add_axes(fr(*CONFIG["layout"]["parameter_icons"]))
profile_ax.set(xlim=(0,3),ylim=(0,1))
profile_ax.axis("off")
# Location: two pins on a ruler.
profile_ax.hlines(.36,.10,.89,color=CONFIG["style"]["location_color"],lw=.65)
profile_ax.vlines(np.linspace(.12,.87,6),.31,.41,color=CONFIG["style"]["location_color"],lw=.5)
profile_ax.vlines([.29,.68],.42,.67,color=CONFIG["style"]["location_color"],lw=.65)
profile_ax.scatter([.29,.68],[.70,.70],s=[16,22],marker="v",
    c=["#C5D5DF",CONFIG["style"]["location_color"]],edgecolors="white",linewidths=.4,zorder=4)
profile_ax.add_patch(FancyArrowPatch((.29,.89),(.68,.89),arrowstyle="->",
    mutation_scale=4.5,lw=.7,color=CONFIG["style"]["location_color"]))
# Scale: a narrow range enclosed by a wider range.
profile_ax.add_patch(FancyBboxPatch((1.11,.40),.78,.28,
    boxstyle="round,pad=0,rounding_size=.07",facecolor="#DCEBE6",
    edgecolor=CONFIG["style"]["scale_color"],lw=.65))
profile_ax.add_patch(FancyBboxPatch((1.30,.44),.40,.20,
    boxstyle="round,pad=0,rounding_size=.05",facecolor=CONFIG["style"]["scale_color"],
    edgecolor="white",lw=.4))
profile_ax.vlines(1.5,.31,.78,color=OUTLINE,lw=.45,linestyles=":")
profile_ax.add_patch(FancyArrowPatch((1.12,.88),(1.88,.88),arrowstyle="|-|",
    mutation_scale=3.5,lw=.7,color=CONFIG["style"]["scale_color"]))
# Shape: asymmetric nested quantile spans, with no density curve.
fan_probabilities=np.column_stack([np.arange(1,5)/10,1-np.arange(1,5)/10])
fan_quantiles=np.exp(mu+sigma/Q*np.log(gammaincinv(k,fan_probabilities)/k))
fan_x=2.10+.78*(fan_quantiles-fan_quantiles.min())/np.ptp(fan_quantiles)
fan_y=np.linspace(.38,.87,4)
profile_ax.hlines(fan_y,fan_x[:,0],fan_x[:,1],color=CONFIG["style"]["shape_color"],lw=1.05)
profile_ax.scatter(fan_x.ravel(),np.repeat(fan_y,2),s=5,
    color=CONFIG["style"]["shape_color"],edgecolors="white",linewidths=.2,zorder=4)
fan_median=np.exp(mu+sigma/Q*np.log(gammaincinv(k,.5)/k))
profile_ax.vlines(2.10+.78*(fan_median-fan_quantiles.min())/np.ptp(fan_quantiles),
    .31,.94,color=OUTLINE,lw=.5,linestyles=":")
parameter_labels=np.frompyfunc(profile_ax.text,3,1)(
    [.5,1.5,2.5],[.02,.02,.02],[r"$\mu$",r"$\sigma$",r"$Q$"])
plt.setp(parameter_labels.tolist(),ha="center",fontsize=9.5)
art.text(1.99,.99,r"$C^2=\frac{1+C_+^2}{\pi}-1$",ha="center",fontsize=10.1)
art.text(1.99,.72,"Total-abundance CV²",ha="center",fontsize=7.2,color=OUTLINE,fontfamily="Inter")

# Paired states carry the core across-gene computation.
px,py,pw,ph=CONFIG["layout"]["population"]
art.add_patch(FancyBboxPatch((px+depth_x,py+depth_y),pw,ph,
    boxstyle="round,pad=.018,rounding_size=.10",facecolor=CONFIG["style"]["board_side"],
    edgecolor=CONFIG["style"]["board_edge"],lw=.55,zorder=.35))
population_board=FancyBboxPatch((px,py),pw,ph,boxstyle="round,pad=.018,rounding_size=.10",
    facecolor=CONFIG["style"]["board_front"],edgecolor=CONFIG["style"]["board_edge"],lw=.6,zorder=.5)
population_board.set_path_effects([pe.SimplePatchShadow(offset=(1.4,-1.5),alpha=.10,
    shadow_rgbFace=CONFIG["style"]["shadow_color"]),pe.Normal()])
art.add_patch(population_board)
art.plot([px+.08,px+pw-.08],[py+ph-.018,py+ph-.018],color="#FFFFFF",lw=1.0,zorder=.6)
art.text(3.50,3.52,r"$L_{g_1}(k),\ldots,L_{g_G}(k)$",ha="center",fontsize=10.6)
connector((3.50,3.39),(3.50,3.22),OUTLINE)
state_ax=fig.add_axes(fr(*CONFIG["layout"]["state_grid"]))
weights=np.asarray(CONFIG["schematic"]["weights"])
weights=(weights+weights.T)/2
sns.heatmap(weights,ax=state_ax,cmap=sns.blend_palette(["#FFFCF2","#E8D6A5","#9DC9C7"],as_cmap=True),
    vmin=0,vmax=weights.max(),cbar=False,xticklabels=False,yticklabels=False,linewidths=2.2,linecolor="white")
columns=np.tile(np.arange(3),3)
rows=np.repeat(np.arange(3),3)
state_vectors=np.column_stack([CONFIG["schematic"]["state_sigma"],
    CONFIG["schematic"]["state_Q"],CONFIG["schematic"]["state_epsilon"],
    CONFIG["schematic"]["state_relative_shift"],CONFIG["schematic"]["state_nonzero"]])
fingerprints=state_vectors/state_vectors.max(axis=0)
fingerprint_y=rows[:,None]+.30+np.arange(5)[None,:]*.13
fingerprint_left=np.broadcast_to(columns[:,None]+.11,(9,5))
fingerprint_right=np.broadcast_to(columns[:,None]+.57,(9,5))
state_ax.hlines(fingerprint_y.ravel(),fingerprint_left.ravel(),
    (fingerprint_left+.30).ravel(),color=MEAN,alpha=.13,lw=1.7,capstyle="round")
state_ax.hlines(fingerprint_y.ravel(),fingerprint_right.ravel(),
    (fingerprint_right+.30).ravel(),color=DV,alpha=.13,lw=1.7,capstyle="round")
state_ax.hlines(fingerprint_y.ravel(),fingerprint_left.ravel(),
    (fingerprint_left+.30*fingerprints[rows]).ravel(),color=MEAN,lw=1.7,capstyle="round")
state_ax.hlines(fingerprint_y.ravel(),fingerprint_right.ravel(),
    (fingerprint_right+.30*fingerprints[columns]).ravel(),color=DV,lw=1.7,capstyle="round")
state_pair_labels=np.char.add(np.char.add((rows+1).astype(str)," | "),(columns+1).astype(str))
state_pair_text=np.frompyfunc(state_ax.text,3,1)(columns+.5,rows+.15,state_pair_labels)
plt.setp(state_pair_text.tolist(),ha="center",va="center",fontsize=5.5,color=OUTLINE,fontfamily="Inter")
state_ax.set(xlim=(0,3),ylim=(3,0))
art.text(3.50,1.91,"Rows: control · cols: case",ha="center",fontsize=7.1,fontfamily="Inter")
prior_ax=fig.add_axes(fr(*CONFIG["layout"]["prior"]))
prior_ax.bar(np.arange(9),weights.ravel(),width=.65,color="#77ABA5",edgecolor="white",linewidth=.35)
prior_ax.axis("off")
art.text(3.50,1.41,r"Shared prior $\widehat V_k$",ha="center",fontsize=7.5,fontfamily="Inter")
art.text(3.50,1.14,r"$w_{gk}\propto L_g(k)\widehat V_k$",ha="center",fontsize=10.0)
art.text(3.50,.76,"Separate mean / CV² states",ha="center",fontsize=7.1,color=OUTLINE,fontfamily="Inter")

# Quantile objects illustrate isolated mean or CV2 changes.
contrast_k=np.array(CONFIG["schematic"]["contrast_shape"])[None,:]
contrast_mean=np.array(CONFIG["schematic"]["contrast_mean"])[None,:]
contrast_scale=contrast_mean/contrast_k
mean_ax=fig.add_axes(fr(*CONFIG["layout"]["mean_icon"]))
cloud=CONFIG["schematic"]["cloud"]
cloud_probabilities=np.linspace(cloud["lower"],cloud["upper"],cloud["points"])[:,None]
cloud_quantiles=gammaincinv(contrast_k[:,:2],cloud_probabilities)*contrast_scale[:,:2]
cloud_offsets=np.where(np.arange(cloud["points"])%2==0,.055,-.055)
mean_ax.hlines([.33,.75],cloud_quantiles.min(axis=0),cloud_quantiles.max(axis=0),
    colors=[MEAN,DV],lw=5,alpha=.12,capstyle="round")
mean_ax.scatter(cloud_quantiles[:,0],.33+cloud_offsets,s=8,color=MEAN,
    edgecolors="white",linewidths=.25,zorder=3)
mean_ax.scatter(cloud_quantiles[:,1],.75+cloud_offsets,s=8,color=DV,
    edgecolors="white",linewidths=.25,zorder=3)
mean_ax.scatter(contrast_mean[0,:2],[.33,.75],s=29,marker="D",c=[MEAN,DV],
    edgecolors="white",linewidths=.7,zorder=4)
mean_ax.vlines(contrast_mean[0,:2],[.42,.84],1.01,colors=[MEAN,DV],lw=.6,linestyles=":")
mean_ax.annotate("",xy=(contrast_mean[0,1],1.05),xytext=(contrast_mean[0,0],1.05),
    arrowprops={"arrowstyle":"->","color":OUTLINE,"lw":.8,"mutation_scale":5})
mean_ax.set(xlim=(0,3.6),ylim=(0,1.18))
mean_ax.axis("off")
art.text(5.16,3.67,"Mean",ha="center",fontsize=8.3,fontweight="semibold",color=MEAN)
art.text(5.16,2.49,r"$\Delta_M=\log_2\frac{E(X_1)}{E(X_0)}$",ha="center",fontsize=10.2)
cv_ax=fig.add_axes(fr(*CONFIG["layout"]["cv_icon"]))
band_quantiles=gammaincinv(contrast_k[:,2:4],
    np.asarray(CONFIG["schematic"]["band_probabilities"])[:,None])*contrast_scale[:,2:4]
cv_ax.barh([.75,.33],band_quantiles[3]-band_quantiles[0],height=.24,
    left=band_quantiles[0],color=[MEAN,DV],edgecolor=[MEAN,DV],linewidth=.7,alpha=.24)
cv_ax.barh([.75,.33],band_quantiles[2]-band_quantiles[1],height=.24,
    left=band_quantiles[1],color=[MEAN,DV],edgecolor="none",alpha=.65)
cv_ax.scatter(contrast_mean[0,2:4],[.75,.33],s=24,marker="D",c=[MEAN,DV],
    edgecolors="white",linewidths=.7,zorder=4)
cv_ax.axvline(CONFIG["schematic"]["contrast_mean"][2],color=OUTLINE,lw=.55,ls=":",zorder=0)
cv_ax.set(xlim=(0,4.5),ylim=(0,1.18))
cv_ax.axis("off")
art.text(5.16,2.21,"Relative variability",ha="center",fontsize=8.3,fontweight="semibold",color=DV)
art.text(5.16,1.02,r"$\Delta_V=\log_2(C_1^2/C_0^2)$",ha="center",fontsize=10.1)

# A folded vector document encodes all four two-channel decisions.
dx,dy,dw,dh=CONFIG["layout"]["document"]
document=Path([(dx,dy),(dx+dw,dy),(dx+dw,dy+dh-.18),(dx+dw-.18,dy+dh),
    (dx,dy+dh),(dx,dy)],
    [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY])
art.add_patch(FancyBboxPatch((dx+depth_x*1.4,dy+depth_y*1.4),dw,dh-.04,
    boxstyle="round,pad=0,rounding_size=.025",facecolor="#C4D5CA",edgecolor="#A8BEB2",lw=.5,zorder=1.3))
art.add_patch(FancyBboxPatch((dx+depth_x*.7,dy+depth_y*.7),dw,dh-.02,
    boxstyle="round,pad=0,rounding_size=.025",facecolor=CONFIG["style"]["paper_back"],
    edgecolor="#B5C9BE",lw=.5,zorder=1.5))
document_front=PathPatch(document,facecolor=CONFIG["style"]["paper_front"],edgecolor="#86A5AE",lw=.85,zorder=2)
document_front.set_path_effects([pe.SimplePatchShadow(offset=(1.4,-1.4),alpha=.11,
    shadow_rgbFace=CONFIG["style"]["shadow_color"]),pe.Normal()])
art.add_patch(document_front)
art.add_patch(FancyBboxPatch((dx+.045,dy+.05),dw-.09,dh-.10,
    boxstyle="round,pad=0,rounding_size=.035",facecolor="#F7FBFA",edgecolor="none",zorder=2.1))
art.add_patch(Rectangle((dx+.055,dy+dh-.055),dw-.30,.018,
    facecolor=MEAN,edgecolor="none",alpha=.75,zorder=3))
art.add_patch(Polygon([(dx+dw-.18,dy+dh),(dx+dw-.18,dy+dh-.18),(dx+dw,dy+dh-.18)],
    facecolor="#F4E7C4",edgecolor="#86A5AE",lw=.65,zorder=3))
art.text(6.61,3.37,"Per gene",ha="center",fontsize=7.8,fontweight="semibold")
art.text(6.29,3.11,"M",ha="center",fontsize=7.0,color=MEAN,fontweight="semibold")
art.text(6.47,3.11,"V",ha="center",fontsize=7.0,color=DV,fontweight="semibold")
art.text(6.80,3.11,"Call",ha="center",fontsize=7.0)
call_y=np.array([2.84,2.47,2.10,1.73])
art.hlines(call_y-.17,6.20,7.01,color="#E8EFF0",lw=.6,zorder=2)
art.scatter(np.repeat(6.29,4),call_y,s=16,c=[MEAN,"white",MEAN,"white"],
    edgecolors=MEAN,linewidths=.65,zorder=4)
art.scatter(np.repeat(6.47,4),call_y,s=16,c=["white",DV,DV,"white"],
    edgecolors=DV,linewidths=.65,zorder=4)
call_text=np.frompyfunc(art.text,3,1)(np.repeat(6.61,4),call_y,["DE only","DV only","Both","Neither"])
plt.setp(call_text.tolist(),fontsize=7.1,va="center",ha="left")
art.text(6.61,1.39,"Effects + intervals",ha="center",fontsize=6.9,fontfamily="Inter")
art.text(6.61,.81,"BH per channel",ha="center",fontsize=7.4)
art.text(6.61,.57,r"$\alpha=0.05$",ha="center",fontsize=10.0)

# Arrows connect the named input, model and output boundaries.
links=CONFIG["layout"]["links"]
tuple(map(connector,np.asarray(links)[:,0,:],np.asarray(links)[:,1,:],
    [OUTLINE,OUTLINE,MEAN,DV,MEAN,DV]))
art.text(5.94,3.35,r"$q_M$",ha="center",fontsize=9.5,color=MEAN)
art.text(5.94,1.92,r"$q_V$",ha="center",fontsize=9.5,color=DV)
art.set(xlim=(0,W),ylim=(0,H),aspect="equal")
fig.savefig(F / "Figure_1.svg")
fig.savefig(F / "Figure_1.pdf")
fig.savefig(F / "Figure_1.png",dpi=450)
plt.close(fig)
plt.rcParams.update(previous_plot_style)
