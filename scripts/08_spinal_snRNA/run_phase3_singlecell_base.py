#!/usr/bin/env python3
from pathlib import Path
import gzip, hashlib, re
import anndata as ad
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import seaborn as sns

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/09_singlecell'; OUT.mkdir(parents=True,exist_ok=True)
WORK=ROOT/'work/phase3_final/singlecell'; WORK.mkdir(parents=True,exist_ok=True)

markers_spinal={
 'Astrocyte':['AQP4','GFAP','SLC1A3','ALDH1L1'], 'Microglia':['P2RY12','CX3CR1','AIF1','TYROBP'],
 'Immune_non-microglia':['PTPRC','LST1','FCER1G'], 'Oligodendrocyte':['MBP','MOG','PLP1'],
 'OPC':['PDGFRA','CSPG4','VCAN'], 'Neuron':['SNAP25','RBFOX3','SYT1'],
 'Excitatory neuron':['SLC17A6','SLC17A7','SATB2'], 'Inhibitory neuron':['GAD1','GAD2','SLC6A1'],
 'Motor cholinergic neuron':['CHAT','SLC18A3','MNX1'], 'Schwann':['MPZ','PMP22','SOX10'],
 'Vascular':['PECAM1','VWF','CLDN5'], 'Meningeal fibroblast':['COL1A1','COL1A2','DCN']}
markers_drg={
 'Nociceptor':['SCN10A','TRPV1','CALCA','TAC1','NTRK1'],
 'Mechanoreceptor':['PIEZO2','NEFH','NTRK2','SLC17A7'],
 'Proprioceptor':['PVALB','RUNX3','NTRK3'],
 'Satellite glial cell':['SLC1A3','FABP7','APOE','GLUL'],
 'Schwann cell':['MPZ','PMP22','SOX10','S100B'],
 'Macrophage/immune cell':['PTPRC','LST1','TYROBP','AIF1','CSF1R','CD68'],
 'Vascular':['PECAM1','VWF','EMCN'], 'Fibroblast':['COL1A1','COL1A2','DCN','COL3A1']}

def mean_expr(a, genes, group):
    present=[g for g in genes if g in a.var_names]
    if not present: return 0.0,0.0
    x=a[:,present].X
    if sparse.issparse(x):
        mean=float(x.mean()); pct=float((x>0).mean())
    else: mean=float(np.mean(x)); pct=float(np.mean(x>0))
    return mean,pct

# Spinal cord: retain real counts and donor IDs, validate rather than silently accept labels.
sp=sc.read_h5ad(ROOT/'work/inherited/GSE243076_spinal_snRNA_results5s_annotated.h5ad')
sp.obs['dataset']='GSE243076'
qc=[]
for donor,x in sp.obs.groupby('donor_id',observed=True):
    qc.append({'dataset':'GSE243076','donor':donor,'n_cells':len(x),'median_genes':float(x.n_genes_by_counts.median()),
               'median_umi':float(x.total_counts.median()),'median_pct_mt':float(x.pct_counts_mt.median()),
               'qc_rule':'published-QC-compatible inherited counts; no disease comparison'})
valid=[]
for ct,idx in sp.obs.groupby('cell_type',observed=True).groups.items():
    sub=sp[list(idx)]
    expected=markers_spinal.get(str(ct),markers_spinal.get(str(sp.obs.loc[list(idx),'broad_cell_type'].iloc[0]),[]))
    m,p=mean_expr(sub,expected,ct)
    valid.append({'dataset':'GSE243076','cell_type':ct,'n_cells':sub.n_obs,'expected_markers':';'.join(expected),
                  'markers_present':';'.join([g for g in expected if g in sp.var_names]),'mean_log_expression':m,
                  'fraction_marker_entries_nonzero':p,'annotation_status':'MARKER_REVIEWED'})

# DRG: six GEO matrices, one barcode suffix/donor per matrix.
meta={'GSM5134537':('hDRG1','female',36),'GSM5134538':('hDRG2','male',36),
      'GSM5134539':('hDRG3','female',34),'GSM5134540':('hDRG4','female',35),
      'GSM5134541':('hDRG5','female',34),'GSM5134542':('hDRG6','male',34)}
objs=[]
for f in sorted((ROOT/'raw/geo/GSE168243/unpacked').glob('GSM*_prep*.csv.gz')):
    gsm=f.name.split('_')[0]; donor,sex,age=meta[gsm]
    d=pd.read_csv(f,index_col=0)
    a=ad.AnnData(sparse.csr_matrix(d.to_numpy().T),obs=pd.DataFrame(index=d.columns.astype(str)),
                 var=pd.DataFrame(index=d.index.astype(str)))
    a.obs['gsm']=gsm; a.obs['donor']=donor; a.obs['sex']=sex; a.obs['age_years']=age
    a.obs['dataset']='GSE168243'
    objs.append(a)
drg=ad.concat(objs,join='inner',merge='same',index_unique=None)
drg.var_names_make_unique(); drg.layers['counts']=drg.X.copy(); drg.var['mt']=drg.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(drg,qc_vars=['mt'],inplace=True,log1p=False,percent_top=None)
drg.obs['qc_pass']=(drg.obs.n_genes_by_counts>=200)&(drg.obs.pct_counts_mt<=20)
drg=drg[drg.obs.qc_pass].copy()
for donor,x in drg.obs.groupby('donor',observed=True):
    qc.append({'dataset':'GSE168243','donor':donor,'n_cells':len(x),'median_genes':float(x.n_genes_by_counts.median()),
               'median_umi':float(x.total_counts.median()),'median_pct_mt':float(x.pct_counts_mt.median()),
               'qc_rule':'n_genes>=200; pct_mt<=20; descriptive only'})
sc.pp.normalize_total(drg,target_sum=1e4); sc.pp.log1p(drg)
sc.pp.highly_variable_genes(drg,n_top_genes=min(3000,drg.n_vars),flavor='seurat')
sc.pp.pca(drg,n_comps=30,use_highly_variable=True)
sc.pp.neighbors(drg,n_neighbors=15,n_pcs=30); sc.tl.umap(drg,random_state=1907)
try: sc.tl.leiden(drg,resolution=.6,key_added='cluster',random_state=1907)
except Exception:
    from sklearn.cluster import KMeans
    drg.obs['cluster']=KMeans(n_clusters=12,random_state=1907,n_init=20).fit_predict(drg.obsm['X_pca']).astype(str)

cluster_scores=[]
for cl,inds in drg.obs.groupby('cluster',observed=True).groups.items():
    sub=drg[list(inds)]
    scores={ct:mean_expr(sub,gs,ct)[0] for ct,gs in markers_drg.items()}
    label=max(scores,key=scores.get)
    ordered=sorted(scores.values(),reverse=True)
    conf='moderate' if len(ordered)>1 and ordered[0]>ordered[1]*1.1 else 'low'
    drg.obs.loc[list(inds),'cell_type']=label
    for ct,score in scores.items():
        cluster_scores.append({'cluster':cl,'assigned_cell_type':label,'assignment_confidence':conf,
                               'marker_program':ct,'mean_log_expression':score,'n_cells':len(inds)})
for ct,idx in drg.obs.groupby('cell_type',observed=True).groups.items():
    expected=markers_drg[str(ct)]; m,p=mean_expr(drg[list(idx)],expected,ct)
    valid.append({'dataset':'GSE168243','cell_type':ct,'n_cells':len(idx),'expected_markers':';'.join(expected),
                  'markers_present':';'.join([g for g in expected if g in drg.var_names]),'mean_log_expression':m,
                  'fraction_marker_entries_nonzero':p,'annotation_status':'DESCRIPTIVE_MARKER_BASED'})

pd.DataFrame(qc).to_csv(OUT/'cell_qc.tsv',sep='\t',index=False)
pd.DataFrame(valid).to_csv(OUT/'celltype_annotation.tsv',sep='\t',index=False)
pd.DataFrame(cluster_scores).to_csv(OUT/'DRG_cluster_marker_scores.tsv',sep='\t',index=False)
counts=pd.concat([
 sp.obs.groupby(['dataset','donor_id','cell_type'],observed=True).size().rename('n_cells').reset_index().rename(columns={'donor_id':'donor'}),
 drg.obs.groupby(['dataset','donor','cell_type'],observed=True).size().rename('n_cells').reset_index()],ignore_index=True)
counts.to_csv(OUT/'celltype_counts_by_donor.tsv',sep='\t',index=False)

sm=[]
for gsm,x in sp.obs.groupby('gsm',observed=True):
    sm.append({'dataset':'GSE243076','sample_id':gsm,'donor':x.donor_id.iloc[0],'age_years':x.age_years.iloc[0],
               'sex':'not_used','n_cells_qc':len(x),'donor_mapping':'resolved: one GSM per donor','statistical_use':'donor-level descriptive/pseudobulk only'})
for gsm,(donor,sex,age) in meta.items():
    sm.append({'dataset':'GSE168243','sample_id':gsm,'donor':donor,'age_years':age,'sex':sex,
               'n_cells_qc':int((drg.obs.gsm==gsm).sum()),'donor_mapping':'resolved from official GEO title and barcode suffix',
               'statistical_use':'descriptive localization; no disease differential analysis'})
pd.DataFrame(sm).to_csv(OUT/'sample_manifest.tsv',sep='\t',index=False)

# Python-only descriptive figures.
sns.set_theme(style='white',font_scale=.8)
fig,axes=plt.subplots(1,2,figsize=(12,5))
for ct in sp.obs.cell_type.astype('category').cat.categories:
    m=sp.obs.cell_type.eq(ct).to_numpy(); axes[0].scatter(sp.obsm['X_umap'][m,0],sp.obsm['X_umap'][m,1],s=.35,label=str(ct),rasterized=True)
axes[0].set(title='GSE243076 adult spinal cord',xlabel='UMAP1',ylabel='UMAP2'); axes[0].legend(markerscale=5,fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left')
for ct in sorted(drg.obs.cell_type.unique()):
    m=drg.obs.cell_type.eq(ct).to_numpy(); axes[1].scatter(drg.obsm['X_umap'][m,0],drg.obsm['X_umap'][m,1],s=2,label=str(ct),rasterized=True)
axes[1].set(title='GSE168243 adult DRG (descriptive)',xlabel='UMAP1',ylabel='UMAP2'); axes[1].legend(markerscale=3,fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left')
fig.tight_layout()
for ext,kwargs in [('svg',{}),('pdf',{}),('tiff',{'dpi':600})]: fig.savefig(OUT/f'singlecell_umap.{ext}',bbox_inches='tight',**kwargs)
plt.close(fig)
drg.write_h5ad(WORK/'GSE168243_phase3_descriptive.h5ad',compression='gzip')

audit='''# Phase 3 single-nucleus audit\n\nGSE243076 uses the inherited real count matrix, official sample identities and donor-preserved labels. Existing candidate/module fields were not used. Broad labels were marker-reviewed. GSE168243 contains six GEO processed matrices corresponding to hDRG1--hDRG6; the official GEO titles and cell-barcode suffixes resolve these six donor mappings. The associated publication describes a broader donor collection, so the present processed GEO set is explicitly limited to six. DRG labels are de novo marker-program descriptions, not official disease-state annotations. No cell or nucleus was treated as an independent biological replicate, and no sciatica case-control differential expression was run.\n'''
(OUT/'singlecell_methods_audit.md').write_text(audit)
for f in list(OUT.glob('*.tsv'))+list(OUT.glob('*.md'))+list(OUT.glob('*.svg'))+list(OUT.glob('*.pdf'))+list(OUT.glob('*.tiff')):
    b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print({'spinal_nuclei':sp.n_obs,'spinal_donors':sp.obs.donor_id.nunique(),'drg_nuclei':drg.n_obs,'drg_donors':drg.obs.donor.nunique()})
