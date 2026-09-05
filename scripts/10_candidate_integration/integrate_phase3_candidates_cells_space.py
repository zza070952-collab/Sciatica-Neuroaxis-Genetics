#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/09_singlecell'
C=pd.read_csv(ROOT/'results/phase3_final/07_genetic_candidates/provisional_genetic_candidates.tsv',sep='\t')
genes=C.gene.astype(str).tolist() if 'gene' in C else C.gene_symbol.astype(str).tolist()

def candidate_expr(h5ad,donor_col,type_col,dataset):
    a=ad.read_h5ad(h5ad)
    if 'paper_qc_pass' in a.obs: a=a[a.obs.paper_qc_pass.astype(bool) & (~a.obs.predicted_doublet.astype(bool))].copy()
    elif 'qc_pass' in a.obs: a=a[a.obs.qc_pass.astype(bool)].copy()
    idx=pd.Index(a.var_names.astype(str)).get_indexer(genes); present=np.where(idx>=0)[0]
    X=a.layers['counts'] if 'counts' in a.layers else a.X; lib=np.asarray(X.sum(1)).ravel()
    rows=[]
    for j in present:
        v=X[:,idx[j]]; v=np.asarray(v.toarray()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
        log=np.log1p(v/np.maximum(lib,1)*1e4)
        tmp=pd.DataFrame({'donor':a.obs[donor_col].astype(str).to_numpy(),'cell_type':a.obs[type_col].astype(str).to_numpy(),'log_expr':log,'detected':v>0})
        for (donor,ct),d in tmp.groupby(['donor','cell_type']): rows.append(dict(dataset=dataset,gene=genes[j],donor=donor,cell_type=ct,n_cells=len(d),mean_log_expression=d.log_expr.mean(),fraction_expressing=d.detected.mean()))
    return pd.DataFrame(rows)

spinal=candidate_expr(ROOT/'work/inherited/GSE243076_spinal_snRNA_results5s_annotated.h5ad','donor_id','cell_type','GSE243076')
drg=candidate_expr(ROOT/'work/phase3_final/singlecell/GSE168243_phase3_descriptive.h5ad','donor','cell_type','GSE168243')
spinal.to_csv(OUT/'candidate_expression_by_celltype.tsv',sep='\t',index=False)
drg.to_csv(OUT/'DRG_candidate_expression.tsv',sep='\t',index=False)

# Trait spatial z versus reference-cell similarity, summarized by section and donor.
M=pd.read_csv(OUT/'spinal_spatial_celltype_mapping.tsv',sep='\t')
S=pd.read_csv(ROOT/'results/phase3_final/08_gsmap/spot_level_results.tsv.gz',sep='\t')
Z=S[S.trait.eq('SCIATICA')][['sample','donor','spot','z']]
X=M.merge(Z,on=['sample','donor','spot'],how='inner'); sim=[c for c in M.columns if c.startswith('similarity_')]
assoc=[]
for (sample,donor),d in X.groupby(['sample','donor']):
    for col in sim:
        rho,p=spearmanr(d[col],d.z)
        assoc.append(dict(sample=sample,donor=donor,cell_type=col.replace('similarity_','').replace('_',' '),spearman_rho=rho,spot_level_p=p,n_spots=len(d)))
A=pd.DataFrame(assoc); A.to_csv(OUT/'gsmap_celltype_association.tsv',sep='\t',index=False)

# Candidate expression localization in each section; spots are descriptive spatial observations.
align=[]
for (sample,trait),z in S.groupby(['sample','trait']):
    a=ad.read_h5ad(ROOT/f'work/phase3_final/gsmap_input_pca/{sample}.h5ad')
    layer=a.layers['count'] if 'count' in a.layers else a.X; lib=np.asarray(layer.sum(1)).ravel(); smap=pd.Series(np.arange(a.n_obs),index=a.obs_names.astype(str))
    zz=z.set_index('spot').loc[a.obs_names.astype(str)]; cutoff=np.quantile(zz.z,.95)
    for gene in genes:
        if gene not in a.var_names: continue
        j=a.var_names.get_loc(gene); v=layer[:,j]; v=np.asarray(v.toarray()).ravel() if sp.issparse(v) else np.asarray(v).ravel(); e=np.log1p(v/np.maximum(lib,1)*1e4)
        rho,p=spearmanr(e,zz.z); high=zz.z.to_numpy()>=cutoff
        align.append(dict(sample=sample,donor=z.donor.iloc[0],trait=trait,gene=gene,n_spots=len(e),fraction_expressing=(v>0).mean(),spearman_rho=rho,spearman_p=p,
          mean_logexpr_top5=e[high].mean(),mean_logexpr_other=e[~high].mean(),top5_minus_other=e[high].mean()-e[~high].mean()))
AL=pd.DataFrame(align); AL.to_csv(OUT/'candidate_spatial_alignment.tsv',sep='\t',index=False)

for f in [OUT/'candidate_expression_by_celltype.tsv',OUT/'DRG_candidate_expression.tsv',OUT/'gsmap_celltype_association.tsv',OUT/'candidate_spatial_alignment.tsv']:
    b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print({'spinal_rows':len(spinal),'drg_rows':len(drg),'celltype_spatial_rows':len(A),'candidate_spatial_rows':len(AL)})
