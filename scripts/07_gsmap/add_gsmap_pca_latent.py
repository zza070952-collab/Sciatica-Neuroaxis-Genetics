#!/usr/bin/env python3
"""Add a documented PCA latent representation to each independent section."""
from pathlib import Path
import anndata as ad
import scanpy as sc
ROOT=Path('.')
src=ROOT/'work/phase3_final/gsmap_input'; dst=ROOT/'work/phase3_final/gsmap_input_pca'; dst.mkdir(parents=True,exist_ok=True)
for f in sorted(src.glob('GSM*.h5ad')):
    a=ad.read_h5ad(f)
    q=a.copy()
    sc.pp.normalize_total(q,target_sum=1e4); sc.pp.log1p(q)
    sc.pp.highly_variable_genes(q,n_top_genes=min(3000,q.n_vars),flavor='seurat')
    sc.pp.pca(q,n_comps=50,use_highly_variable=True,zero_center=False,random_state=1907)
    a.obsm['X_pca_phase3']=q.obsm['X_pca'].copy()
    a.uns['X_pca_phase3_provenance']={'normalization':'library-size 1e4 + log1p','HVG':'seurat 3000','components':50,'random_state':1907,
                                      'purpose':'official gsMap --latent_representation option; raw counts retained in layer count'}
    a.write_h5ad(dst/f.name,compression='gzip')
    print(f.name,a.shape,a.obsm['X_pca_phase3'].shape,flush=True)
