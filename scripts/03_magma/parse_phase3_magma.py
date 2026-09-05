#!/usr/bin/env python3
from pathlib import Path
import hashlib
import math
import re
import numpy as np
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/03_magma'
RAW=OUT/'raw'

loc_path=next((ROOT/'work/software/phase3_magma/NCBI37.3').glob('*.gene.loc'))
loc=pd.read_csv(loc_path,sep=r'\s+',header=None,
                names=['gene','gene_chr','gene_start','gene_stop','strand','gene_symbol'],dtype={'gene':str})
loc['gene_chr']=pd.to_numeric(loc['gene_chr'],errors='coerce').astype('Int64')

def bh(p):
    p=np.asarray(p,float); n=len(p); order=np.argsort(p); q=np.empty(n); prev=1.0
    for rank,idx in reversed(list(enumerate(order,1))):
        prev=min(prev,p[idx]*n/rank); q[idx]=prev
    return q

genes=[]
for trait in ['SCIATICA','LDH']:
    for analysis,suffix in [('primary',''),('no_MHC','_noMHC')]:
        f=RAW/f'{trait}_10kb{suffix}.genes.out'
        d=pd.read_csv(f,sep=r'\s+',comment='#')
        d.columns=[x.lower() for x in d.columns]
        d['gene']=d['gene'].astype(str)
        d=d.merge(loc,on='gene',how='left',validate='many_to_one')
        d['mhc_flag']=(d.gene_chr.eq(6)&d.gene_start.lt(34_000_000)&d.gene_stop.gt(25_000_000))
        d.insert(0,'trait',trait); d.insert(1,'analysis',analysis)
        genes.append(d)
g=pd.concat(genes,ignore_index=True)
g['fdr']=g.groupby(['trait','analysis'])['p'].transform(lambda x: bh(x.to_numpy()))
g['bonferroni_threshold']=g.groupby(['trait','analysis'])['p'].transform(lambda x: .05/len(x))
g.to_csv(OUT/'gene_results.tsv',sep='\t',index=False)
g.loc[(g.analysis=='primary') & ((g.fdr<.05)|(g.p<g.bonferroni_threshold))].to_csv(
    OUT/'gene_results_FDR.tsv',sep='\t',index=False)

def read_gsa(path,trait,collection,analysis):
    lines=path.read_text().splitlines()
    start=next(i for i,x in enumerate(lines) if re.match(r'^VARIABLE\s+',x))
    data=[]
    for line in lines[start:]:
        if not line.strip() or line.startswith('#'): continue
        fields=line.split()
        if fields[0]=='VARIABLE': header=[x.lower() for x in fields]; continue
        if len(fields)!=len(header): continue
        data.append(fields)
    d=pd.DataFrame(data,columns=header)
    for c in d.columns:
        if c not in ['variable','type','full_name']: d[c]=pd.to_numeric(d[c],errors='coerce')
    d.insert(0,'trait',trait); d.insert(1,'analysis',analysis); d.insert(2,'collection',collection)
    return d

sets=[]
for trait in ['SCIATICA','LDH']:
    for label in ['C2','GO']:
        sets.append(read_gsa(RAW/f'{trait}_{label}.gsa.out',trait,label,'primary'))
        nf=RAW/f'{trait}_{label}_noMHC.gsa.out'
        if nf.exists(): sets.append(read_gsa(nf,trait,label,'no_MHC'))
s=pd.concat(sets,ignore_index=True)
s['fdr']=s.groupby(['trait','analysis','collection'])['p'].transform(lambda x: bh(x.fillna(1).to_numpy()))
s.to_csv(OUT/'pathway_results.tsv',sep='\t',index=False)
s.loc[(s.fdr<.05)].to_csv(OUT/'pathway_results_FDR.tsv',sep='\t',index=False)

prim=g[g.analysis.eq('primary')].copy()
cmp=prim.pivot(index=['gene','gene_symbol','gene_chr','gene_start','gene_stop','mhc_flag'],columns='trait',values=['p','fdr','zstat']).reset_index()
cmp.columns=['_'.join([str(y) for y in x if str(y)]) if isinstance(x,tuple) else x for x in cmp.columns]
cmp.to_csv(OUT/'sciatica_LDH_gene_comparison.tsv',sep='\t',index=False)

report=f'''# Phase 3 MAGMA analysis audit

- Reference: 1000 Genomes European LD, GRCh37, MAGMA v1.10.
- Variant alignment: rsID intersection; FinnGen ALT remained the effect allele in upstream QC.
- Primary gene window: 10 kb upstream/downstream.
- Tested genes: SCIATICA {len(prim[prim.trait.eq('SCIATICA')]):,}; LDH {len(prim[prim.trait.eq('LDH')]):,}.
- Bonferroni thresholds: SCIATICA {0.05/len(prim[prim.trait.eq('SCIATICA')]):.3g}; LDH {0.05/len(prim[prim.trait.eq('LDH')]):.3g}.
- Pathways: MSigDB 2026.1 C2 canonical (including Reactome/KEGG where present) and C5 GO.
- Sensitivity: chr6:25--34 Mb removed before gene analysis.
- LDH is a structural-reference phenotype and cannot nominate a SCIATICA candidate by itself.
'''
(OUT/'magma_methods_audit.md').write_text(report)
for f in [OUT/'gene_results.tsv',OUT/'gene_results_FDR.tsv',OUT/'pathway_results.tsv',OUT/'pathway_results_FDR.tsv',OUT/'sciatica_LDH_gene_comparison.tsv',OUT/'magma_methods_audit.md']:
    b=f.read_bytes(); (Path(str(f)+'.md5')).write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); (Path(str(f)+'.sha256')).write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print({'genes':len(g),'fdr_genes':int(((g.analysis=='primary')&(g.fdr<.05)).sum()),'pathways':len(s),'fdr_pathways':int((s.fdr<.05).sum())})
