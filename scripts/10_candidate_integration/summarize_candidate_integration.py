#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/10_candidate_freeze'; OUT.mkdir(parents=True,exist_ok=True)
p=pd.read_csv(ROOT/'results/phase3_final/07_genetic_candidates/provisional_genetic_candidates.tsv',sep='\t')
s=pd.read_csv(ROOT/'results/phase3_final/09_singlecell/candidate_expression_by_celltype.tsv',sep='\t')
d=pd.read_csv(ROOT/'results/phase3_final/09_singlecell/DRG_candidate_expression.tsv',sep='\t')
a=pd.read_csv(ROOT/'results/phase3_final/09_singlecell/candidate_spatial_alignment.tsv',sep='\t')
a=a[a.trait.eq('SCIATICA')]
j=pd.read_csv(ROOT/'results/phase3_final/05_twas/conditional_joint_results.tsv',sep='\t')
j=j[(j.trait.eq('SCIATICA'))&(j.joint_status.eq('INCLUDED'))].copy()
j['ensembl_gene']=j.ID.astype(str).str.split('.').str[0]

def best_cell(frame,gene):
    x=frame[frame.gene.eq(gene)].groupby('cell_type').agg(
        donors=('donor','nunique'),fraction=('fraction_expressing','mean'),mean_expr=('mean_log_expression','mean'),cells=('n_cells','sum')).reset_index()
    if not len(x): return ('NOT_DETECTED',0,0.0,0.0,0)
    r=x.sort_values(['donors','fraction','mean_expr'],ascending=False).iloc[0]
    return (r.cell_type,int(r.donors),r.fraction,r.mean_expr,int(r.cells))

rows=[]
for r in p.itertuples(index=False):
    sx=best_cell(s,r.gene); dx=best_cell(d,r.gene)
    ax=a[a.gene.eq(r.gene)]
    donor=ax.groupby('donor').spearman_rho.median()
    rows.append({
        'gene':r.gene,'tier':r.tier,'MAGMA_FDR':r.SCIATICA_MAGMA_FDR,
        'TWAS_min_FDR':r.SCIATICA_TWAS_min_FDR,'TWAS_robust_tissues':r.SCIATICA_TWAS_robust_tissues,
        'TWAS_joint_included':bool((j.ensembl_gene==str(r.ensembl_gene).split('.')[0]).any()),
        'ABF_max_PP4':r.coloc_max_PP4,'SuSiE_max_PP4':r.susie_max_PP4,
        'LDH_role':r.LDH_role,'LDH_MAGMA_FDR':r.LDH_MAGMA_FDR,
        'spatial_median_rho':ax.spearman_rho.median(),'spatial_positive_sections':int((ax.spearman_rho>0).sum()),
        'spatial_valid_sections':int(ax.spearman_rho.notna().sum()),'spatial_positive_donors':int((donor>0).sum()),
        'spinal_top_cell':sx[0],'spinal_donors':sx[1],'spinal_fraction':sx[2],
        'DRG_top_cell':dx[0],'DRG_donors':dx[1],'DRG_fraction':dx[2]
    })
o=pd.DataFrame(rows).sort_values(['tier','SuSiE_max_PP4','ABF_max_PP4','MAGMA_FDR'],ascending=[True,False,False,True],na_position='last')
path=OUT/'candidate_integration_summary.tsv'; o.to_csv(path,sep='\t',index=False)
b=path.read_bytes(); Path(str(path)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+path.name+'\n'); Path(str(path)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+path.name+'\n')
print(o.to_string(index=False))
