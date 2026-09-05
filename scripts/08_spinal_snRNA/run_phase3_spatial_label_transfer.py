#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/09_singlecell'; OUT.mkdir(parents=True,exist_ok=True)
REF=ROOT/'work/inherited/GSE243076_spinal_snRNA_results5s_annotated.h5ad'
SP=ROOT/'work/phase3_final/gsmap_input_pca'
qc=pd.read_csv(ROOT/'results/phase3_final/08_gsmap/section_qc.tsv',sep='\t')

a=ad.read_h5ad(REF)
keep=a.obs['paper_qc_pass'].astype(bool).to_numpy() & (~a.obs['predicted_doublet'].astype(bool).to_numpy())
keep &= ~a.obs['cell_type'].isin(['Ambiguous broad identity','Motor/cholinergic neuron']).to_numpy()
a=a[keep].copy(); X=a.layers['counts'] if 'counts' in a.layers else a.X
if not sp.issparse(X): X=sp.csr_matrix(X)
X=X.tocsr(); labels=a.obs.cell_type.astype(str).to_numpy(); donors=a.obs.donor_id.astype(str).to_numpy()
types=sorted(np.unique(labels)); genes=pd.Index(a.var_names.astype(str))

def aggregate(indices, group_labels, groups):
    indicator=sp.csr_matrix((np.ones(len(indices)),(group_labels,indices)),shape=(len(groups),X.shape[0]))
    return indicator @ X

type_idx=np.array([types.index(x) for x in labels])
pb=aggregate(np.arange(X.shape[0]),type_idx,types)
pb=np.asarray(pb.todense(),dtype=np.float64)
pb=np.log1p(pb/np.maximum(pb.sum(1,keepdims=True),1)*1e6)

# Select de novo Phase 3 markers by pseudobulk specificity, not inherited candidate/module scores.
bad=genes.str.startswith(('MT-','RPL','RPS'))
marker_by={}
for i,t in enumerate(types):
    other=np.max(np.delete(pb,i,axis=0),axis=0)
    score=pb[i]-other; score[bad]=-np.inf
    ix=np.argsort(score)[-100:]; ix=ix[np.isfinite(score[ix]) & (score[ix]>0)]
    marker_by[t]=genes[ix].tolist()
marker_union=pd.Index(sorted(set(sum(marker_by.values(),[]))))
gix=genes.get_indexer(marker_union)

def zrows(mat):
    mat=np.asarray(mat,dtype=np.float32)
    mu=mat.mean(1,keepdims=True); sd=mat.std(1,keepdims=True); sd[sd<1e-6]=1
    return (mat-mu)/sd

refz=zrows(pb[:,gix])

# Donor-level leave-one-donor-out validation on donor x cell-type pseudobulks.
cv=[]
for held in sorted(np.unique(donors)):
    train=donors!=held; test=donors==held
    train_pb=[]
    for t in types:
        v=np.asarray(X[train & (labels==t)].sum(0)).ravel()
        train_pb.append(np.log1p(v/max(v.sum(),1)*1e6))
    train_pb=np.vstack(train_pb)
    fold_markers=[]
    for j,t in enumerate(types):
        other=np.max(np.delete(train_pb,j,axis=0),axis=0)
        score=train_pb[j]-other; score[bad]=-np.inf
        ix=np.argsort(score)[-100:]; ix=ix[np.isfinite(score[ix]) & (score[ix]>0)]
        fold_markers.extend(genes[ix].tolist())
    fold_union=pd.Index(sorted(set(fold_markers))); fold_gix=genes.get_indexer(fold_union)
    trainz=zrows(train_pb[:,fold_gix])
    for true in types:
        if not np.any(test & (labels==true)): continue
        v=np.asarray(X[test & (labels==true)].sum(0)).ravel(); v=np.log1p(v/max(v.sum(),1)*1e6)
        vz=zrows(v[None,fold_gix]); cor=(vz @ trainz.T / len(fold_gix)).ravel(); pred=types[int(np.argmax(cor))]
        cv.append(dict(held_out_donor=held,true_cell_type=true,predicted_cell_type=pred,correct=pred==true,
                       top_correlation=float(cor.max()),margin=float(np.sort(cor)[-1]-np.sort(cor)[-2])))
cv=pd.DataFrame(cv); cv.to_csv(OUT/'spinal_label_transfer_LODO_validation.tsv',sep='\t',index=False)

maps=[]
for r in qc.itertuples(index=False):
    s=ad.read_h5ad(SP/f'{r.sample}.h5ad')
    common=marker_union.intersection(pd.Index(s.var_names.astype(str)))
    if len(common)<200: raise RuntimeError(f'too few marker genes for {r.sample}: {len(common)}')
    ridx=marker_union.get_indexer(common); sidx=pd.Index(s.var_names.astype(str)).get_indexer(common)
    sx=s.layers['count'][:,sidx] if 'count' in s.layers else s.X[:,sidx]
    if sp.issparse(sx): lib=np.asarray(sx.sum(1)).ravel(); dense=sx.toarray().astype(np.float32)
    else: dense=np.asarray(sx,dtype=np.float32); lib=dense.sum(1)
    dense=np.log1p(dense/np.maximum(lib[:,None],1)*1e4)
    cor=zrows(dense) @ refz[:,ridx].T / len(common)
    # Softmax scores are similarity weights, not estimated cell fractions.
    e=np.exp((cor-cor.max(1,keepdims=True))/.10); soft=e/e.sum(1,keepdims=True); order=np.argsort(cor,axis=1)
    pred=np.array(types)[order[:,-1]]
    d=pd.DataFrame({'sample':r.sample,'donor':r.donor,'spot':s.obs_names.astype(str),'predicted_cell_type':pred,
                    'top_correlation':cor[np.arange(len(pred)),order[:,-1]],
                    'correlation_margin':cor[np.arange(len(pred)),order[:,-1]]-cor[np.arange(len(pred)),order[:,-2]],
                    'n_marker_genes_used':len(common)})
    for j,t in enumerate(types): d['similarity_'+t.replace(' ','_').replace('/','_')]=soft[:,j]
    maps.append(d)
mapping=pd.concat(maps,ignore_index=True)
mapping.to_csv(OUT/'spinal_spatial_celltype_mapping.tsv',sep='\t',index=False)

markers=[]
for t,g in marker_by.items():
    for rank,gene in enumerate(reversed(g),1): markers.append((t,rank,gene))
pd.DataFrame(markers,columns=['cell_type','specificity_rank','gene']).to_csv(OUT/'spinal_reference_markers.tsv',sep='\t',index=False)
summary=pd.DataFrame([{'metric':'LODO_donor_celltype_accuracy','value':cv.correct.mean(),'n_units':len(cv)},
                      {'metric':'spatial_spots_mapped','value':len(mapping),'n_units':len(mapping)},
                      {'metric':'reference_cell_types','value':len(types),'n_units':len(types)},
                      {'metric':'marker_union_genes','value':len(marker_union),'n_units':len(marker_union)}])
summary.to_csv(OUT/'spinal_label_transfer_summary.tsv',sep='\t',index=False)

audit=f'''# Spinal reference-to-Visium label-transfer audit

- Reference: GSE243076 adult human spinal-cord snRNA-seq; paper-QC nuclei retained and predicted doublets removed.
- Candidate/module columns present in the inherited audited object were not used. Phase 3 reference markers were recalculated from count-layer donor-aware pseudobulk profiles.
- Validation unit: donor x cell type. Leave-one-donor-out reference rebuilding yielded accuracy {cv.correct.mean():.3f} across {len(cv)} held-out pseudobulk units.
- Mapping: correlation label transfer over {len(marker_union)} de novo cell-type-specific genes shared with each Visium section. Softmax-transformed correlations are similarity weights, not cell fractions and not disease-state estimates.
- Sections were mapped separately. Donors remain the biological replicates; no case-control differential expression was performed.
'''
(OUT/'spinal_label_transfer_audit.md').write_text(audit)
for f in [OUT/'spinal_label_transfer_LODO_validation.tsv',OUT/'spinal_spatial_celltype_mapping.tsv',OUT/'spinal_reference_markers.tsv',OUT/'spinal_label_transfer_summary.tsv',OUT/'spinal_label_transfer_audit.md']:
    b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print(summary.to_string(index=False))
