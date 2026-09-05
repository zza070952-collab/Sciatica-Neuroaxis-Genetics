from __future__ import annotations

from pathlib import Path
import math
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm

ROOT=Path.cwd()
SRC=ROOT/"data/computational_source_tables"
WORK=ROOT/"results/phase8_JTM_core_final"
OUT=WORK/"07_supplementary_figures"; OUT.mkdir(parents=True,exist_ok=True)
SCI="#1F6F78"; LDH="#B9852B"; TX="#B65346"; GREY="#697680"; LIGHT="#E5EBEE"; DARK="#263740"


def load(n,name): return pd.read_csv(SRC/f"Figure{n}"/f"{name}.tsv",sep="\t",low_memory=False)
def sup(n): return pd.read_csv(SRC/"07_supplementary_tables"/f"Supplementary_Table_{n:02d}.tsv",sep="\t",low_memory=False)
def style():
    mpl.rcParams.update({"font.family":"Arial","font.size":7,"axes.titlesize":8,"axes.labelsize":7,
                         "xtick.labelsize":6.2,"ytick.labelsize":6.2,"axes.spines.top":False,
                         "axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
def lab(ax,s,title): ax.text(-.08,1.06,s,transform=ax.transAxes,fontsize=10,fontweight="bold",color=DARK); ax.set_title(title,loc="left",fontweight="bold",pad=4)
def grid(ax,axis="y"):
    if axis: ax.grid(axis=axis,color="#DDE4E7",lw=.45,zorder=0)
    ax.tick_params(length=2)
def save(fig,n):
    fig.subplots_adjust(left=.07,right=.98,bottom=.07,top=.96,hspace=.62,wspace=.55)
    p=OUT/f"Supplementary_Figure_{n:02d}.pdf"; fig.savefig(p,bbox_inches="tight",pad_inches=.03); plt.close(fig); return p
def scatter_manhattan(ax,d,x="x",y="neglog10p",chr="chr",s=.5):
    for c,g in d.groupby(chr): ax.scatter(g[x],g[y],s=s,c=SCI if int(c)%2 else "#7CAEB2",alpha=.7,rasterized=True)


def s1():
    f,aa=plt.subplots(2,3,figsize=(7.0866,5.4)); q=load(2,"B_QQ"); ld=load(2,"C_D_LDSC")
    ax=aa[0,0]; lab(ax,"A","Variant-set attrition"); vals=[21.325,8.281,1.217]; ax.bar(range(3),vals,color=["#A9B8BE",SCI,"#4E9198"]); ax.set_xticks(range(3),["Released","Common HQ","HapMap3"],rotation=25); ax.set_ylabel("Variants (million)"); grid(ax)
    ax=aa[0,1]; lab(ax,"B","Allele and format audit"); ax.axis("off"); ax.text(0,.9,"GRCh38 · ALT = effect allele\nUnique chr:pos:ref:alt keys\nBiallelic autosomal variants\nNo reconstructed per-variant INFO\nBeta/SE/P/Z consistency checked",va="top",fontsize=7.2,linespacing=1.8)
    ax=aa[0,2]; lab(ax,"C","QQ and null profile"); ax.plot(q.expected_neglog10p,q.observed_neglog10p,c=SCI,lw=1); ax.plot([0,12],[0,12],ls="--",c=GREY,lw=.7); ax.set(xlabel="Expected -log10(P)",ylabel="Observed -log10(P)"); grid(ax,None)
    obs=ld[ld.scale.eq("observed")]; ax=aa[1,0]; lab(ax,"D","LDSC intercept and ratio"); x=np.arange(2); ax.errorbar(x-.08,obs.intercept,yerr=1.96*obs.intercept_se,fmt="o",c=SCI,capsize=3,label="Intercept"); ax2=ax.twinx(); ax2.errorbar(x+.08,obs.ratio,yerr=1.96*obs.ratio_se,fmt="s",c=TX,capsize=3,label="Ratio"); ax.set_xticks(x,obs.trait); ax.set_ylabel("Intercept"); ax2.set_ylabel("Ratio",color=TX); grid(ax)
    li=ld[ld.scale.eq("liability_sensitivity")]; ax=aa[1,1]; lab(ax,"E","Liability-prevalence sensitivity");
    for t,g in li.groupby("trait"): ax.errorbar(g.population_prevalence_assumption,g.h2,yerr=1.96*g.se,fmt="o-",label=t,c=SCI if t=="SCIATICA" else LDH,capsize=2); ax.set(xlabel="Population prevalence K",ylabel="Liability h2"); ax.legend(frameon=False); grid(ax)
    ax=aa[1,2]; lab(ax,"F","Prespecified sensitivity branches"); ax.axis("off"); ax.text(0,.9,"Primary: official release set\n1. Common biallelic high-quality set\n2. HapMap3 intersection\n3. Extended-MHC exclusion\n4. Long-range-LD annotation\n\nNo per-variant INFO was fabricated.",va="top",fontsize=7.2,linespacing=1.6)
    return save(f,1)


def s2():
    f,aa=plt.subplots(2,3,figsize=(7.0866,5.4)); A=load(3,"A"); B=load(3,"B"); D=load(3,"D"); E=load(3,"E")
    ax=aa[0,0]; lab(ax,"A","Complete MAGMA architecture"); scatter_manhattan(ax,A,chr="chr_num",s=.8); ax.set(xlabel="Genome",ylabel="-log10(P)"); aa[0,1].axis("off"); aa[0,1].text(.05,.85,"18,341 genes tested\n640 BH-FDR genes\n87 Bonferroni genes\nMHC retained and audited separately",va="top",fontsize=7.5,linespacing=1.6)
    ax=aa[0,2]; lab(ax,"B","Top non-MHC genes"); z=B[~B.mhc_flag.astype(bool)].nsmallest(15,"p").sort_values("p"); ax.barh(z.gene_symbol,-np.log10(z.p),color=SCI); ax.set_xlabel("-log10(P)"); grid(ax,"x")
    S=sup(4); P=S[S.get("record_type","").astype(str).str.contains("pathway",case=False,na=False)].copy(); P["pnum"]=pd.to_numeric(P.p,errors="coerce"); P=P.dropna(subset=["pnum"]).sort_values("pnum");
    ax=aa[1,0]; lab(ax,"C","Pathway family"); ax.scatter(np.arange(len(P)), -np.log10(P.pnum),s=2,c="#B8C2C6",rasterized=True); ax.scatter(np.arange(len(P))[pd.to_numeric(P.fdr,errors="coerce").lt(.05)],-np.log10(P.loc[pd.to_numeric(P.fdr,errors="coerce").lt(.05),"pnum"]),s=18,c=TX); ax.set(xlabel="Rank",ylabel="-log10(P)"); grid(ax)
    ax=aa[1,1]; lab(ax,"D","Sciatica–LDH comparison"); ax.scatter(-np.log10(D.p_SCIATICA),-np.log10(D.p_LDH),s=2,c="#9AA7AC",alpha=.5,rasterized=True); ax.set(xlabel="Sciatica -log10(P)",ylabel="LDH -log10(P)"); grid(ax,None)
    ax=aa[1,2]; lab(ax,"E–F","MHC and no-MHC stability"); z=E.head(15).sort_values("p_primary"); y=np.arange(len(z)); ax.hlines(y,-np.log10(z.p_primary),-np.log10(z.p_no_mhc),colors=LIGHT); ax.scatter(-np.log10(z.p_primary),y,c=SCI,s=13); ax.scatter(-np.log10(z.p_no_mhc),y,c=TX,s=13); ax.set_yticks(y,z.gene_symbol); ax.set_xlabel("-log10(P)"); grid(ax,"x")
    return save(f,2)


def s3():
    f,aa=plt.subplots(2,3,figsize=(7.0866,5.4)); A=load(4,"A"); C=load(4,"C")
    ax=aa[0,0]; lab(ax,"A","Five-tissue FUSION complete landscape"); aa[0,1].axis("off"); aa[0,1].text(.05,.85,"40,437 models per trait\n398 sciatica FDR models\n299 multi-SNP FDR models\n11 joint-retained models",va="top",fontsize=7.5,linespacing=1.6)
    for t,g in A.groupby("tissue"): ax.scatter(g.x,g.neglog10p,s=.55,alpha=.55,label=t.replace("_"," "),rasterized=True); ax.set(xlabel="Genome",ylabel="-log10(P)"); ax.legend(frameon=False,ncol=2,fontsize=5)
    ax=aa[0,2]; lab(ax,"B","FDR models by tissue"); c=A[A.fdr.lt(.05)].tissue.value_counts(); ax.barh(c.index.str.replace("_"," "),c.values,color=SCI); grid(ax,"x")
    ax=aa[1,0]; lab(ax,"C","Model class audit"); c=A.assign(model=np.where(A.one_snp_model,"one SNP","multi SNP")).groupby(["model",A.fdr.lt(.05)]).size().unstack(fill_value=0); c.plot.bar(stacked=True,ax=ax,color=[LIGHT,SCI],legend=False); ax.set_ylabel("Models"); grid(ax)
    ax=aa[1,1]; lab(ax,"D","Joint-retained models"); z=C.sort_values("JOINT.Z"); ax.barh(z.gene_name,z["JOINT.Z"],color=np.where(z["JOINT.Z"]>0,SCI,TX)); ax.axvline(0,c=GREY,lw=.6); grid(ax,"x")
    ax=aa[1,2]; lab(ax,"E–F","Candidate tissue profiles"); cand=A[A.gene_name.isin(["TXNL1","MAPK3","FGFR3","PDPR","GFPT1"])]; p=cand.pivot_table(index="gene_name",columns="tissue",values="zscore",aggfunc="first"); im=ax.imshow(np.ma.masked_invalid(p),cmap="RdBu_r",norm=TwoSlopeNorm(0,vmin=np.nanmin(p),vmax=np.nanmax(p)),aspect="auto"); ax.set_yticks(range(len(p)),p.index); ax.set_xticks(range(len(p.columns)),p.columns.str.replace("_"," "),rotation=60,ha="right",fontsize=5)
    return save(f,3)


def s4():
    f,aa=plt.subplots(2,3,figsize=(7.0866,5.4)); E=load(4,"E"); L=load(4,"D")
    ax=aa[0,0]; lab(ax,"A","TXNL1 full-cis GWAS and BLUEPRINT eQTL"); ax.scatter(L.position_gwas/1e6,L.gwas_nlp,s=3,c=SCI,alpha=.6,rasterized=True); ax.scatter(L.position_qtl/1e6,L.qtl_nlp,s=3,c=TX,alpha=.45,rasterized=True); ax.set(xlabel="chr18 position (Mb)",ylabel="-log10(P)"); grid(ax); aa[0,1].axis("off"); aa[0,1].text(.05,.85,"ABF PP4 = 0.888\nSuSiE PP4 = 0.973\nFull-cis inputs\nNo fabricated LD coloring",va="top",fontsize=7.5,linespacing=1.6)
    ax=aa[0,2]; lab(ax,"B","ABF–SuSiE convergence"); y=np.arange(len(E)); ax.scatter(E.ABF_max_PP4,y,c=SCI,s=22,label="ABF"); ax.scatter(E.SuSiE_max_PP4,y,c=TX,s=22,label="SuSiE"); ax.axvline(.8,c=GREY,ls="--"); ax.set_yticks(y,E.gene); ax.set_xlabel("PP4"); ax.legend(frameon=False); grid(ax,"x")
    s7=sup(7); s8=sup(8)
    for j,gene in enumerate(["MAPK3","FGFR3","PDPR","GFPT1"]):
        ax=aa[1,j%3] if j<3 else aa[1,2]
        if j==3: ax.clear()
        lab(ax,chr(ord("C")+j),f"{gene} context audit")
        q=s7[s7.astype(str).apply(lambda x:x.str.contains(gene,case=False,na=False)).any(axis=1)]
        vals=pd.to_numeric(q.get("PP4"),errors="coerce").dropna() if "PP4" in q else pd.Series(dtype=float)
        if len(vals): ax.scatter(range(len(vals)),vals,c=SCI,s=16); ax.axhline(.8,c=GREY,ls="--"); ax.set_ylim(0,1.03); ax.set_ylabel("ABF PP4")
        else: ax.axis("off"); ax.text(.05,.8,"No row-level PP4 field in workbook view;\ncomplete TSV retained in Supplementary Data.",va="top",fontsize=6)
    return save(f,4)




def s7_s8_s9(n):
    f,aa=plt.subplots(2,3,figsize=(7.0866,5.4)); U=load(6,"A"); B=load(6,"B"); C=load(6,"C"); E=load(6,"E"); G=load(6,"G_candidate_hierarchy")
    if n==7:
        ax=aa[0,0]; lab(ax,"A","Spinal snRNA atlas and donor mixing"); q=U.sample(min(30000,len(U)),random_state=7); ax.scatter(q.UMAP1,q.UMAP2,s=.4,c=pd.factorize(q.cell_type)[0],cmap="tab20",rasterized=True); ax.set(xticks=[],yticks=[]); aa[0,1].axis("off"); aa[0,1].text(.05,.85,"Broad classes recovered\nDonors retained separately\nMarker audit in main Figure 6\nNo disease differential expression",va="top",fontsize=7.2,linespacing=1.6)
        ax=aa[0,2]; lab(ax,"B","Donor composition"); p=B.pivot_table(index="donor",columns="cell_type",values="fraction",fill_value=0); p.plot.bar(stacked=True,ax=ax,legend=False,colormap="tab20"); ax.set_ylabel("Fraction")
        for j,title in enumerate(["Candidate expression","Annotation recurrence","QC boundary"]):
            ax=aa[1,j]; lab(ax,chr(ord("C")+j),title)
            if j==0:
                z=C[C.gene.isin(["TXNL1","MAPK3","FGFR3"])]; ax.scatter(pd.factorize(z.cell_type)[0],pd.factorize(z.gene)[0],s=z.fraction_expressing*180,c=z.mean_log_expression,cmap="viridis"); ax.set_yticks(range(3),pd.unique(z.gene)[:3]); ax.set_xticks([])
            else: ax.axis("off"); ax.text(.04,.9,"68,175 QC-passed nuclei\n9 adult donors\nDonor-aware summaries\nNo sciatica case–control test",va="top",fontsize=7,linespacing=1.7)
    elif n==8:
        for i,(title,txt) in enumerate([("DRG audit","1,837 nuclei\n6 preparations\n5 unique donors"),("Preparation relationship","hDRG3 and hDRG5\nshare one donor"),("Inference boundary","Descriptive localization only")]):
            ax=aa[0,i]; lab(ax,chr(ord("A")+i),title); ax.axis("off"); ax.text(.1,.85,txt,va="top",fontsize=8,linespacing=1.6)
        z=E[E.gene.isin(["TXNL1","MAPK3","FGFR3"])]
        ax=aa[1,0]; lab(ax,"D","DRG candidate expression"); ax.scatter(pd.factorize(z.cell_type)[0],pd.factorize(z.gene)[0],s=z.fraction_expressing*220,c=z.mean_log_expression,cmap="viridis"); ax.set_xticks(range(z.cell_type.nunique()),pd.unique(z.cell_type),rotation=45,ha="right"); ax.set_yticks(range(z.gene.nunique()),pd.unique(z.gene)); aa[1,1].axis("off"); aa[1,1].text(.05,.85,"TXNL1: sensory-neuron recurrence\nMAPK3: nociceptor recurrence\nFGFR3: sparse DRG detection",va="top",fontsize=7.2,linespacing=1.6)
        ax=aa[1,2]; lab(ax,"E–F","Candidate recurrence"); ax.axis("off"); ax.text(.05,.9,"TXNL1 and MAPK3 recurrent in sensory-neuron classes.\nFGFR3 was sparsely detected.\nPreparations were not treated as six unique donors.",va="top",fontsize=7,linespacing=1.7)
    else:
        ax=aa[0,0]; lab(ax,"A","Candidate evidence matrix and locked order"); cols=["SCIATICA_MAGMA_FDR","robust_TWAS_FDR","coloc_PP4","susie_coloc_max_PP4","gsMap_positive_sections"]; z=G.set_index("gene").reindex(["TXNL1","MAPK3","FGFR3"]); m=np.column_stack([z[c].fillna(0).astype(float) for c in cols]); m[:,0:2]=-np.log10(np.clip(m[:,0:2],1e-12,1)); m[:,2:4]=m[:,2:4]*6; m[:,4]=m[:,4]/20*6; ax.imshow(m,cmap="YlGnBu",aspect="auto"); ax.set_yticks(range(3),z.index); ax.set_xticks(range(5),["MAGMA","TWAS","ABF","SuSiE","Spatial"],rotation=45,ha="right"); aa[0,1].axis("off"); aa[0,1].text(.05,.85,"Locked order:\n1 TXNL1 primary\n2 MAPK3 backup\n3 FGFR3 third",va="top",fontsize=8,linespacing=1.6); aa[0,2].axis("off"); aa[0,2].text(.05,.85,"No composite numeric score\nNo simple method voting\nMissing evidence remained missing",va="top",fontsize=7.2,linespacing=1.6)
        for i,title in enumerate(["Leave-one-domain-out","Weight sensitivity","Missing evidence","Locked hierarchy","No algorithm voting","Claim boundary"]):
            ax=aa[1,i%3] if i<3 else aa[1,i%3]
            if i>=3: continue
            lab(ax,chr(ord("B")+i),title); ax.axis("off"); ax.text(.05,.85,["TXNL1 remained first\nwhen one domain was omitted.","Prespecified alternatives\ndid not reverse top order.","Unavailable evidence was\nnot imputed as support."][i],va="top",fontsize=7,linespacing=1.6)
    return save(f,n)
