#!/usr/bin/env python3
from pathlib import Path
import re, sqlite3, hashlib
import numpy as np
import pandas as pd

ROOT=Path('.')
MAG=ROOT/'results/phase3_final/03_magma/gene_results.tsv'
TW=ROOT/'results/phase3_final/05_twas/twas_all_results.tsv'
OUT=ROOT/'results/phase3_final/07_genetic_candidates'; OUT.mkdir(parents=True,exist_ok=True)
CO=ROOT/'results/phase3_final/06_coloc'; CO.mkdir(parents=True,exist_ok=True)

m=pd.read_csv(MAG,sep='\t'); m=m[(m.trait.eq('SCIATICA'))&(m.analysis.eq('primary'))].copy()
t=pd.read_csv(TW,sep='\t'); t=t[t.trait.eq('SCIATICA')].copy()
t['ensembl_gene']=t.gene.astype(str).str.split('.').str[0]

# Position and stable gene-ID inventory from the five predeclared PredictDB models.
dbdir=ROOT/'work/software/phase3_twas/predictdb_core/eqtl/mashr'
dbs=sorted(dbdir.glob('*.db'))
rows=[]
for db in dbs:
    tissue=db.stem.replace('mashr_','')
    c=sqlite3.connect(db)
    e=pd.read_sql_query('select gene, genename from extra',c)
    w=pd.read_sql_query('select gene, varID from weights',c)
    c.close()
    x=w.varID.str.extract(r'^chr(?P<chr>[^_]+)_(?P<pos>\d+)_')
    w=pd.concat([w,x],axis=1); w['pos']=pd.to_numeric(w.pos,errors='coerce')
    q=w.groupby('gene').agg(chr=('chr',lambda z:z.mode().iloc[0] if len(z.mode()) else np.nan),
                            center_b38=('pos','median'),model_min_b38=('pos','min'),model_max_b38=('pos','max')).reset_index()
    q=q.merge(e,on='gene',how='left'); q['tissue']=tissue
    q['ensembl_gene']=q.gene.str.split('.').str[0]
    rows.append(q)
dbmap=pd.concat(rows,ignore_index=True)

# GENCODE gene records provide fallback symbol-to-Ensembl mapping for MAGMA-only genes.
gtf=ROOT/'work/phase3_final/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'
gtf_rows=[]
with gtf.open(errors='replace') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f=line.rstrip().split('\t')
        if len(f)<9 or f[2]!='gene': continue
        attrs=dict(re.findall(r'(\S+) "([^"]+)"',f[8]))
        if 'gene_name' in attrs and 'gene_id' in attrs:
            gtf_rows.append((attrs['gene_name'],attrs['gene_id'].split('.')[0],f[0].replace('chr',''),int(f[3]),int(f[4])))
gtfmap=pd.DataFrame(gtf_rows,columns=['gene_symbol','ensembl_gene','gtf_chr','gtf_start_b37','gtf_stop_b37']).drop_duplicates('gene_symbol')

twagg=t.groupby(['gene_name','ensembl_gene']).agg(
    twas_min_p=('pvalue','min'),twas_min_fdr=('fdr','min'),
    twas_fdr_tissues=('fdr',lambda x:int((x<.05).sum())),
    twas_robust_tissues=('n_snps_used',lambda x:int(((x>=2)&(t.loc[x.index,'fdr']<.05)).sum())),
    twas_max_snps=('n_snps_used','max'),twas_tissues=('tissue',lambda x:';'.join(sorted(set(x[t.loc[x.index,'fdr']<.05]))))
).reset_index().rename(columns={'gene_name':'gene_symbol'})
ms=m[['gene','gene_symbol','gene_chr','gene_start','gene_stop','p','fdr','zstat','mhc_flag']].rename(
    columns={'p':'magma_p','fdr':'magma_fdr','zstat':'magma_z'})

pool=pd.merge(ms,twagg,on='gene_symbol',how='outer')
pool['magma_fdr_hit']=pool.magma_fdr.lt(.05)
pool['twas_regulatory_hit']=pool.twas_robust_tissues.fillna(0).ge(1)
pool['twas_multitissue_hit']=pool.twas_robust_tissues.fillna(0).ge(2)
pool=pool[pool.magma_fdr_hit | pool.twas_multitissue_hit].copy()
pool['priority_class']=np.select([
    pool.magma_fdr_hit & pool.twas_regulatory_hit,
    pool.magma_fdr_hit,
    pool.twas_multitissue_hit],['MAGMA_FDR_plus_multisnp_TWAS','MAGMA_FDR','multitissue_multisnp_TWAS'],default='not_eligible')
rank={'MAGMA_FDR_plus_multisnp_TWAS':0,'MAGMA_FDR':1,'multitissue_multisnp_TWAS':2}
pool['_rank']=pool.priority_class.map(rank)
pool=pool.sort_values(['_rank','magma_fdr','twas_min_fdr'],na_position='last').drop_duplicates('gene_symbol')

# Prefer exact GRCh38 model positions; use broad lift37 fallback only for lookup.
exact=dbmap.sort_values(['ensembl_gene','model_min_b38']).groupby('ensembl_gene').agg(
    qtl_chr=('chr','first'),qtl_center=('center_b38','median'),qtl_model_min=('model_min_b38','min'),qtl_model_max=('model_max_b38','max'),
    model_tissues=('tissue',lambda x:';'.join(sorted(set(x))))) .reset_index()
pool=pool.merge(gtfmap,on='gene_symbol',how='left',suffixes=('','_gtf'))
pool['ensembl_gene']=pool.ensembl_gene.combine_first(pool.ensembl_gene_gtf)
pool=pool.merge(exact,on='ensembl_gene',how='left')
gene_chr_fallback=pool.gene_chr.astype('Int64').astype(str).where(pool.gene_chr.notna())
pool['qtl_chr']=pool.qtl_chr.combine_first(pool.gtf_chr).combine_first(gene_chr_fallback)
pool['qtl_center']=pool.qtl_center.combine_first(((pool.gtf_start_b37+pool.gtf_stop_b37)/2)).combine_first(((pool.gene_start+pool.gene_stop)/2))
pool['coordinate_source']=np.where(pool.qtl_model_min.notna(),'PredictDB_GRCh38_weight_positions','GENCODE_lift37_broad_lookup_fallback')
pool['mhc_flag']=pool.mhc_flag.fillna((pool.qtl_chr.astype(str).eq('6'))&(pool.qtl_center.between(25_000_000,34_000_000)))
pool=pool[~pool.mhc_flag].copy().head(20)
pool['coloc_region_start']=(pool.qtl_center-2_000_000).clip(lower=1).round().astype('Int64')
pool['coloc_region_end']=(pool.qtl_center+2_000_000).round().astype('Int64')
pool['selection_status']=np.where(pool.ensembl_gene.notna()&pool.qtl_center.notna(),'TARGETED_COLOC','NOT_AVAILABLE_NO_STABLE_GENE_MAPPING')

cols=['gene_symbol','ensembl_gene','priority_class','magma_p','magma_fdr','magma_z','twas_min_p','twas_min_fdr','twas_fdr_tissues','twas_robust_tissues','twas_max_snps','twas_tissues','qtl_chr','qtl_center','coloc_region_start','coloc_region_end','coordinate_source','mhc_flag','selection_status']
pool[cols].to_csv(CO/'targeted_coloc_candidates_pre_qtl.tsv',sep='\t',index=False)
audit='''# Genetics-only provisional selection audit

This candidate list was generated before reading gsMap or candidate-expression results. Eligibility was based on SCIATICA only: MAGMA FDR < 0.05, or FDR-significant TWAS in at least two predeclared tissues with at least two SNPs used per qualifying model. One-SNP TWAS models did not independently qualify a gene. LDH did not nominate candidates. MHC genes were excluded. The 20-gene cap was applied by a prespecified evidence hierarchy, then by within-method FDR; this is not an algorithm-vote count and is not the final experimental freeze.

Coordinates preferentially use exact GRCh38 PredictDB model variants. A broad +/-2 Mb lookup around GENCODE lift37 coordinates is used only when no PredictDB model coordinate exists; returned QTL rows must still match the exact Ensembl gene ID and exact GWAS/QTL variant alleles.
'''
(OUT/'candidate_selection_audit_pre_coloc.md').write_text(audit)
for f in [CO/'targeted_coloc_candidates_pre_qtl.tsv',OUT/'candidate_selection_audit_pre_coloc.md']:
    b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print(pool[cols].to_string(index=False))
