#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
ROOT=Path('.')
OUT=ROOT/'results/phase3_final/07_genetic_candidates'; OUT.mkdir(parents=True,exist_ok=True)
target=pd.read_csv(ROOT/'results/phase3_final/06_coloc/targeted_coloc_candidates_pre_qtl.tsv',sep='\t')
target=target[['gene_symbol','ensembl_gene','priority_class','qtl_chr','qtl_center','coordinate_source','selection_status']]
m=pd.read_csv(ROOT/'results/phase3_final/03_magma/gene_results.tsv',sep='\t'); t=pd.read_csv(ROOT/'results/phase3_final/05_twas/twas_all_results.tsv',sep='\t'); co=pd.read_csv(ROOT/'results/phase3_final/06_coloc/coloc_all.tsv',sep='\t')
su=pd.read_csv(ROOT/'results/phase3_final/06_coloc/susie_coloc.tsv',sep='\t')

def twagg(trait):
 d=t[t.trait.eq(trait)].copy(); d['robust']=d.fdr.lt(.05)&d.n_snps_used.ge(2)
 return d.groupby('gene_name').agg(**{
  f'{trait}_TWAS_min_p':('pvalue','min'),f'{trait}_TWAS_min_FDR':('fdr','min'),
  f'{trait}_TWAS_FDR_tissues':('fdr',lambda x:int((x<.05).sum())),f'{trait}_TWAS_robust_tissues':('robust','sum'),
  f'{trait}_TWAS_tissues':('tissue',lambda x:';'.join(sorted(set(x[d.loc[x.index,'fdr']<.05]))))}).reset_index().rename(columns={'gene_name':'gene_symbol'})

ms=m[(m.trait.eq('SCIATICA'))&(m.analysis.eq('primary'))][['gene_symbol','gene','gene_chr','gene_start','gene_stop','p','fdr','zstat','mhc_flag']].rename(columns={'gene':'entrez_gene','p':'SCIATICA_MAGMA_p','fdr':'SCIATICA_MAGMA_FDR','zstat':'SCIATICA_MAGMA_z'})
ml=m[(m.trait.eq('LDH'))&(m.analysis.eq('primary'))][['gene_symbol','p','fdr','zstat']].rename(columns={'p':'LDH_MAGMA_p','fdr':'LDH_MAGMA_FDR','zstat':'LDH_MAGMA_z'})
cagg=co.groupby('gene_symbol').agg(coloc_max_PP4=('PP4','max'),coloc_best_context=('context',lambda x:x.iloc[co.loc[x.index,'PP4'].fillna(-1).argmax()] if len(x) else ''),
  coloc_strong_contexts=('context',lambda x:';'.join(sorted(set(x[co.loc[x.index,'PP4']>=.8])))),coloc_completed_contexts=('status',lambda x:int((x=='COMPLETED').sum()))).reset_index()
sagg=su.groupby('gene_symbol').agg(susie_max_PP4=('PP.H4.abf','max'),
  susie_strong_signal_pairs=('PP.H4.abf',lambda x:int((x>=.8).sum())),
  susie_completed_signal_pairs=('status',lambda x:int((x=='COMPLETED').sum())),
  susie_completed_no_signal_pairs=('status',lambda x:int((x=='COMPLETED_NO_SIGNAL_PAIR').sum()))).reset_index()
d=target.rename(columns={'gene_symbol':'gene'}).merge(ms,left_on='gene',right_on='gene_symbol',how='left').drop(columns=['gene_symbol']).merge(ml,left_on='gene',right_on='gene_symbol',how='left').drop(columns=['gene_symbol'])
d=d.merge(twagg('SCIATICA'),left_on='gene',right_on='gene_symbol',how='left').drop(columns=['gene_symbol']).merge(twagg('LDH'),left_on='gene',right_on='gene_symbol',how='left').drop(columns=['gene_symbol']).merge(cagg,left_on='gene',right_on='gene_symbol',how='left').drop(columns=['gene_symbol']).merge(sagg,left_on='gene',right_on='gene_symbol',how='left').drop(columns=['gene_symbol'])
d['MHC_flag']=d.mhc_flag.fillna(False)
d['long_range_LD_flag']=((d.gene_chr.eq(8)&d.gene_start.between(7e6,13e6))|(d.gene_chr.eq(11)&d.gene_start.between(46e6,57e6))|(d.gene_chr.eq(17)&d.gene_start.between(40e6,45e6)))
d['tier']=np.select([
 d.SCIATICA_MAGMA_FDR.lt(.05)&d.coloc_max_PP4.ge(.8),
 d.SCIATICA_MAGMA_FDR.lt(.05)&d.SCIATICA_TWAS_robust_tissues.fillna(0).ge(1)],['G1','G2'],default='G3')
d['tier_reason']=np.select([
 d.tier.eq('G1'),d.tier.eq('G2')],['MAGMA_FDR_plus_strong_full_cis_ABF_coloc','MAGMA_FDR_plus_multisnp_tissue_TWAS'],default='single_or_unconfirmed_strong_genetic_evidence')
d['LDH_role']=np.where(d.LDH_MAGMA_FDR.lt(.05)|d.LDH_TWAS_min_FDR.lt(.05),'structural_reference_concordance','no_strong_LDH_gene_level_support')
order={'G1':0,'G2':1,'G3':2}; d['_rank']=d.tier.map(order); d=d.sort_values(['_rank','coloc_max_PP4','SCIATICA_MAGMA_FDR','SCIATICA_TWAS_min_FDR'],ascending=[True,False,True,True],na_position='last').drop(columns=['_rank'])
d.to_csv(OUT/'provisional_genetic_candidates.tsv',sep='\t',index=False); d.to_csv(OUT/'candidate_evidence_matrix.tsv',sep='\t',index=False)
audit='''# Provisional genetics candidate audit

The genetics-only universe was frozen before candidate expression, gsMap-celltype alignment, or experimental feasibility was inspected. Tier G1 requires SCIATICA MAGMA FDR < 0.05 plus full-cis single-signal ABF coloc PP4 >= 0.80. Tier G2 requires SCIATICA MAGMA FDR < 0.05 plus an FDR-significant tissue TWAS model using at least two SNPs. Tier G3 is a single or not-yet-orthogonally-confirmed strong genetic signal. Multi-signal SuSiE-coloc is retained separately as a robustness diagnostic and does not retroactively redefine the frozen ABF tier. LDH is reported only as structural-reference concordance and cannot nominate a candidate. MHC and predeclared long-range-LD regions are flagged rather than hidden. Ranking is hierarchical by evidence definition, not by counting algorithms.
'''
(OUT/'candidate_selection_audit.md').write_text(audit)
for f in [OUT/'provisional_genetic_candidates.tsv',OUT/'candidate_evidence_matrix.tsv',OUT/'candidate_selection_audit.md']:
 b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n');Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print(d[['gene','tier','SCIATICA_MAGMA_FDR','SCIATICA_TWAS_robust_tissues','coloc_max_PP4']].to_string(index=False))
