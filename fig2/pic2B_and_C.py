#!/usr/bin/env python3
"""Standalone depth-matched Figure 2B/C, with full source data and no P stars."""

from pathlib import Path
import argparse
import hashlib
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/doff_matplotlib_revision")
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REV=Path(__file__).resolve().parents[1]
SAMPLES=tuple(f"{t}_R{r}" for t in ("Root","Stem","Leaf") for r in (1,2,3))
LABELS=tuple(f"{t}{r}" for t in ("Root","Stem","Leaf") for r in (1,2,3))
COLORS=("#32886F",)*3+("#BE6B58",)*3+("#487CA6",)*3


def save(fig,path):
    for ext,kwargs in (("svg",{}),("pdf",{}),("png",{"dpi":250}),("tiff",{"dpi":600})):
        fig.savefig(path.with_suffix('.'+ext),facecolor="white",**kwargs)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=Path,default=REV/"figures/expression_revision")
    parser.add_argument("--overwrite",action="store_true",help="Regenerate these derived panels during visual QA")
    args=parser.parse_args()
    out=args.outdir
    if not args.overwrite and any((out/f"Figure2{p}_depth_matched.svg").exists() for p in "BC"):
        raise ValueError("Do not overwrite completed figure outputs")
    source=REV/"results/20_uniform_count_statistics/aggregate/all_nine_technical_tests_with_global_FDR.tsv"
    frame=pd.read_csv(source,sep="\t")
    frame=frame.loc[frame["testable"]].copy()
    if len(frame)!=2833 or frame.duplicated(["sample","tetrad_id"]).any():
        raise ValueError("Unexpected source population")
    if not np.allclose(frame["plugin_paired_excess"],frame["observed_maximum_fraction"]-frame["plugin_null_median"]):
        raise ValueError("Paired source values do not reconstruct")
    mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
        "font.size":9,"axes.labelsize":10,"xtick.labelsize":9,"ytick.labelsize":9,
        "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.8,
        "axes.grid":False,"svg.fonttype":"none","pdf.fonttype":42})
    out.mkdir(parents=True,exist_ok=True)
    fig_b,ax_b=plt.subplots(figsize=(7.2,3.0))
    fig_c,ax_c=plt.subplots(figsize=(7.2,3.0))
    for fig in (fig_b,fig_c): fig.subplots_adjust(left=.12,right=.99,bottom=.19,top=.89)
    fig_b.subplots_adjust(top=.80)
    summary=[]
    for x,(sample,label,color) in enumerate(zip(SAMPLES,LABELS,COLORS)):
        chosen=frame.loc[frame["sample"]==sample]
        real=chosen["observed_maximum_fraction"].to_numpy()
        null=chosen["plugin_null_median"].to_numpy()
        excess=chosen["plugin_paired_excess"].to_numpy()
        for values,pos,fill in ((null,x-.18,"#B9BDC1"),(real,x+.18,color)):
            violin=ax_b.violinplot(values,positions=[pos],widths=.33,showextrema=False,points=160)
            for body in violin["bodies"]:
                body.set_facecolor(fill); body.set_edgecolor("none"); body.set_alpha(.85)
            low,med,high=np.quantile(values,[.25,.5,.75])
            ax_b.plot([pos,pos],[low,high],color="#222222",lw=1.2,zorder=4)
            ax_b.scatter([pos],[med],s=8,color="#111111",zorder=5)
        ax_b.text(x,1.12,f"n={len(chosen)}",ha="center",va="bottom",fontsize=7.4,transform=ax_b.get_xaxis_transform())
        ax_b.text(x,1.035,f"{np.median(real):.3f}",ha="center",va="bottom",fontsize=7.4,color=color,transform=ax_b.get_xaxis_transform())
        rng=np.random.default_rng(20261103+x)
        ax_c.scatter(x+rng.uniform(-.22,.22,len(excess)),excess,s=5,color=color,alpha=.4,linewidths=0,zorder=2)
        low,med,high=np.quantile(excess,[.25,.5,.75])
        ax_c.plot([x,x],[low,high],color="#111111",lw=2,zorder=4)
        ax_c.plot([x-.13,x+.13],[med,med],color="#111111",lw=1.5,zorder=5)
        ax_c.text(x,1.025,f"{med:.3f}",ha="center",va="bottom",fontsize=7.8,transform=ax_c.get_xaxis_transform())
        summary.append({"sample":sample,"n":len(chosen),"real_median":np.median(real),"null_locus_median":np.median(null),
            "paired_excess_median":med,"paired_excess_q25":low,"paired_excess_q75":high})
    for ax in (ax_b,ax_c):
        ax.set_xlim(-.6,8.6); ax.set_xticks(range(9),LABELS); ax.tick_params(axis="both",direction="out")
        ax.set_facecolor("white")
    ax_b.set_ylim(0,1.01); ax_b.set_yticks([0,.25,.5,.75,1.0],["0","0.25","0.50","0.75","1.00"])
    ax_b.set_ylabel("Maximum homolog fraction")
    values=frame["plugin_paired_excess"].to_numpy()
    low=min(-.05,float(values.min())-.025); high=float(values.max())+.035
    ax_c.set_ylim(low,high); ax_c.axhline(0,color="#555555",lw=.8,ls="--",zorder=1)
    ax_c.set_ylabel("Observed \u2212 null median",fontname="DejaVu Sans")
    ax_b.text(-.10,1.08,"B",transform=ax_b.transAxes,fontsize=14,fontweight="bold")
    ax_c.text(-.10,1.08,"C",transform=ax_c.transAxes,fontsize=14,fontweight="bold")
    save(fig_b,out/"Figure2B_depth_matched")
    save(fig_c,out/"Figure2C_depth_matched")
    winner_source=REV/"results/20_uniform_count_statistics/aggregate/all_tissue_winner_sensitivity.tsv"
    winner=pd.read_csv(winner_source,sep="\t")
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.8))
    fig.subplots_adjust(left=.10,right=.99,bottom=.20,top=.75,wspace=.30)
    for ax,subset in zip(axes,("Primary","Primary_clean")):
        chosen=winner.loc[(winner["subset"]==subset)&(winner["class"]=="3_of_3_same")]
        for i,tissue in enumerate(("Root","Stem","Leaf")):
            row=chosen.loc[(chosen["tissue"]==tissue)&(chosen["method"]=="plugin")].iloc[0]
            predictive=chosen.loc[(chosen["tissue"]==tissue)&(chosen["method"]=="dirichlet_predictive")].iloc[0]
            ax.scatter(i-.17,row["observed_pct"],marker="D",s=25,c="#333333",label="Observed" if i==0 else None)
            ax.scatter(i,row["null_expected_pct"],marker="o",s=26,c="#7F9CAB",label="Plug-in null" if i==0 else None)
            ax.scatter(i+.17,predictive["null_expected_pct"],marker="s",s=24,c="#BF875B",label="Predictive null" if i==0 else None)
            ax.text(i,1.02,f"n={int(row['n_loci'])}",ha="center",va="bottom",fontsize=7.5,transform=ax.get_xaxis_transform())
        ax.set_title("Primary" if subset=="Primary" else "Initial Primary-clean",fontsize=9,pad=22)
        ax.set_ylim(0,35);ax.set_xlim(-.6,2.6);ax.set_xticks(range(3),("Root","Stem","Leaf"))
        ax.set_yticks((0,10,20,30))
    axes[0].set_ylabel("Same largest copy in 3 libraries (%)",fontsize=9)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.55,1.02),ncol=3,frameon=False,fontsize=8)
    save(fig,out/"FigureS7_winner_consistency")
    winner.to_csv(out/"FigureS7_source_data.tsv",sep="\t",index=False)
    frame[["sample","tetrad_id","observed_maximum_fraction","plugin_null_median","plugin_paired_excess"]].to_csv(
        out/"Figure2BC_source_data.tsv",sep="\t",index=False)
    pd.DataFrame(summary).to_csv(out/"Figure2BC_summary.tsv",sep="\t",index=False)
    (out/"figure_source_manifest.json").write_text(json.dumps({"source":str(source),
        "sha256":hashlib.sha256(source.read_bytes()).hexdigest(),"threshold":10,"Primary_fixed":534,
        "B_gray_distribution":"Across-locus distribution of locus-specific depth-matched null medians, not pooled null draws",
        "C_points":"All testable locus/library records; black median and IQR; not biological confidence intervals",
        "winner_source":str(winner_source),"winner_source_sha256":hashlib.sha256(winner_source.read_bytes()).hexdigest(),
        "new_statistics_computed":False},indent=2)+"\n")
    print(f"DONE: standalone B/C panels and matched source data: {out}")


if __name__=="__main__": main()
