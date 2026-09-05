#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/10_candidate_freeze'; OUT.mkdir(parents=True,exist_ok=True)
P=pd.read_csv(ROOT/'results/phase3_final/07_genetic_candidates/provisional_genetic_candidates.tsv',sep='\t')
S=pd.read_csv(ROOT/'results/phase3_final/09_singlecell/candidate_expression_by_celltype.tsv',sep='\t')
D=pd.read_csv(ROOT/'results/phase3_final/09_singlecell/DRG_candidate_expression.tsv',sep='\t')
A=pd.read_csv(ROOT/'results/phase3_final/09_singlecell/candidate_spatial_alignment.tsv',sep='\t')
A=A[A.trait.eq('SCIATICA')]
TW=pd.read_csv(ROOT/'results/phase3_final/05_twas/twas_all_results.tsv',sep='\t')
J=pd.read_csv(ROOT/'results/phase3_final/05_twas/conditional_joint_results.tsv',sep='\t')
joint_ids=set(J[(J.trait.eq('SCIATICA'))&(J.joint_status.eq('INCLUDED'))].ID.astype(str).str.split('.').str[0])

def top_cell(df,gene):
    x=df[df.gene.eq(gene)].groupby('cell_type').agg(mean_fraction=('fraction_expressing','mean'),mean_logexpr=('mean_log_expression','mean'),n_donors=('donor','nunique'),n_cells=('n_cells','sum')).reset_index()
    if not len(x): return ('NOT_DETECTED',np.nan,np.nan,0,0)
    # Prefer donor-replicated localization over a higher fraction from a
    # single-donor cell class (notably the DRG Schwann-cell record).
    x=x.sort_values(['n_donors','mean_fraction','mean_logexpr'],ascending=False).iloc[0]
    return (x.cell_type,x.mean_fraction,x.mean_logexpr,int(x.n_donors),int(x.n_cells))

rows=[]
for r in P.itertuples(index=False):
    sct,sfrac,sexpr,sdon,scells=top_cell(S,r.gene)
    dct,dfrac,dexpr,ddon,dcells=top_cell(D,r.gene)
    x=A[A.gene.eq(r.gene)]; dm=x.groupby('donor').spearman_rho.median()
    robust=TW[(TW.trait.eq('SCIATICA'))&(TW.gene_name.eq(r.gene))&(TW.fdr.lt(.05))&(TW.n_snps_used.ge(2))].sort_values('pvalue')
    rb=robust.iloc[0] if len(robust) else None
    rows.append(dict(gene=r.gene,ensembl_gene=r.ensembl_gene,tier=r.tier,
      primary_method='MAGMA+FUSION_multi-SNP_TWAS+full-cis_coloc',SCIATICA_MAGMA_p=r.SCIATICA_MAGMA_p,SCIATICA_MAGMA_FDR=r.SCIATICA_MAGMA_FDR,
      robust_TWAS_tissue=(rb.tissue if rb is not None else 'NONE'),robust_TWAS_z=(rb.zscore if rb is not None else np.nan),robust_TWAS_p=(rb.pvalue if rb is not None else np.nan),robust_TWAS_FDR=(rb.fdr if rb is not None else np.nan),
      FUSION_joint_included=str(r.ensembl_gene).split('.')[0] in joint_ids,
      coloc_PP4=r.coloc_max_PP4,coloc_context=r.coloc_best_context,susie_coloc_max_PP4=r.susie_max_PP4,
      LDH_MAGMA_p=r.LDH_MAGMA_p,LDH_MAGMA_FDR=r.LDH_MAGMA_FDR,LDH_role=r.LDH_role,MHC_flag=r.MHC_flag,long_range_LD_flag=r.long_range_LD_flag,
      gsMap_region='relative high-z lumbar-spinal neighborhoods; no official horn label',gsMap_candidate_median_rho=x.spearman_rho.median(),gsMap_positive_sections=int((x.spearman_rho>0).sum()),gsMap_valid_sections=int(x.spearman_rho.notna().sum()),gsMap_positive_donors=int((dm>0).sum()),
      snRNA_spinal_celltype=sct,snRNA_spinal_fraction=sfrac,snRNA_spinal_donors=sdon,DRG_celltype=dct,DRG_fraction=dfrac,DRG_donors=ddon,
      pain_pseudobulk='NOT_PERFORMED_NO_CASE_CONTROL_DESIGN',direct_SNP_signal='NOT_ESTIMATED_CORE_PIPELINE',required_assay='RNA_KNOCKDOWN_AND_DOWNSTREAM_RESPONSE'))
R=pd.DataFrame(rows)

# Frozen after genetics, gsMap and single-cell completion. These ranks are an
# auditable multi-axis decision, not a numerical algorithm vote.
rank={'TXNL1':1,'MAPK3':2,'FGFR3':3}
R['final_rank']=R.gene.map(rank)
R['final_status']=np.select([R.final_rank.eq(1),R.final_rank.eq(2),R.final_rank.eq(3)],['PRIMARY','BACKUP','THIRD_PREEXPERIMENT'],default='NOT_ADVANCED')
R['cell_model_feasibility']=np.select([R.gene.eq('TXNL1'),R.gene.eq('MAPK3'),R.gene.eq('FGFR3')],
 ['HIGH_human_iPSC_sensory_neuron_with_delivery_optimization','HIGH_human_iPSC_sensory_neuron_with_MAPK1_compensation_control','HIGH_human_stem_cell_derived_spinal_astrocyte'],default='NOT_ASSESSED_BEYOND_TOP3')
R=R.sort_values(['final_rank','tier','SCIATICA_MAGMA_FDR'],na_position='last')
R.to_csv(OUT/'candidate_registry_final.tsv',sep='\t',index=False)

print(R[["gene", "final_rank", "final_status"]].to_string(index=False))
