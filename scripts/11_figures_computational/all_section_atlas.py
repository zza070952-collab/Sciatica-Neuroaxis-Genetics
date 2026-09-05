#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/08_gsmap'
FIG=OUT/'spatial_figures'; FIG.mkdir(parents=True,exist_ok=True)
spots=pd.read_csv(OUT/'spot_level_results.tsv.gz',sep='\t')
sec=pd.read_csv(OUT/'section_level_results.tsv',sep='\t')
don=pd.read_csv(OUT/'donor_level_results.tsv',sep='\t')
loo=pd.read_csv(OUT/'leave_one_donor_out.tsv',sep='\t')
qc=pd.read_csv(OUT/'section_qc.tsv',sep='\t')
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8,'axes.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})

def save(fig,stem,w,h):
    fig.set_size_inches(w,h)
    fig.savefig(FIG/f'{stem}.svg',bbox_inches='tight')
    fig.savefig(FIG/f'{stem}.pdf',bbox_inches='tight')
    fig.savefig(FIG/f'{stem}.tiff',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)

def section_grid(trait):
    d=spots[spots.trait.eq(trait)]; lim=float(np.nanquantile(np.abs(d.z),.98)); lim=max(lim,1)
    fig,axs=plt.subplots(5,4,sharex=False,sharey=False)
    for ax,(sample,x) in zip(axs.ravel(),d.groupby('sample',sort=True)):
        ax.scatter(x.spatial_x,x.spatial_y,c=x.z,s=2.5,cmap='RdBu_r',vmin=-lim,vmax=lim,rasterized=True)
        ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        s=sec[(sec['sample'].eq(sample))&(sec.trait.eq(trait))].iloc[0]
        ax.set_title(f'{sample} | {s.donor}\nACAT P={s.global_p_cauchy:.2g}; FDR spots={s.n_fdr_0_05}')
    sm=mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-lim,lim),cmap='RdBu_r')
    fig.colorbar(sm,ax=axs.ravel().tolist(),fraction=.012,pad=.01,label='gsMap z')
    fig.suptitle(f'{trait}: all 20 normal adult lumbar spinal cord sections',y=.995,fontsize=11)
    save(fig,f'{trait}_all_20_sections',10,12)


section_grid('SCIATICA')
