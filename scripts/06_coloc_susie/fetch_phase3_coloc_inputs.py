#!/usr/bin/env python3
from pathlib import Path
import subprocess, gzip, hashlib, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/06_coloc'; W=ROOT/'work/phase3_final/coloc/inputs'
W.mkdir(parents=True,exist_ok=True)
CAND=OUT/'targeted_coloc_candidates_pre_qtl.tsv'
TAB='tabix'
GW=str(ROOT/'work/phase3_final/coloc/SCIATICA.coloc_eaf.tsv.bgz')
qcols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split()
gcols='rsid chr position ref effect_allele z beta se p EAF'.split()
contexts=[
 ('spinal_cord','QTD000201',126,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000201/QTD000201.all.tsv.gz'),
 ('fibroblast','QTD000216',483,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000216/QTD000216.all.tsv.gz'),
 ('skeletal_muscle','QTD000281',702,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000281/QTD000281.all.tsv.gz'),
 ('tibial_nerve','QTD000286',532,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000286/QTD000286.all.tsv.gz'),
 ('whole_blood','QTD000356',670,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz'),
 ('monocyte_BLUEPRINT','QTD000021',191,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000021/QTD000021.all.tsv.gz'),
 ('macrophage_naive','QTD000001',84,'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000001/QTD000001/QTD000001.all.tsv.gz')]

def run_tabix(path,region,retries=3):
    errs=[]
    for k in range(retries):
        p=subprocess.run([TAB,path,region],text=True,capture_output=True)
        if p.returncode==0: return p.stdout,''
        errs.append(p.stderr.strip())
        time.sleep(1)
    return '', ' | '.join(errs)

c=pd.read_csv(CAND,sep='\t'); c=c[c.selection_status.eq('TARGETED_COLOC')].copy()
c['qtl_chr']=c.qtl_chr.astype(str).str.replace('.0','',regex=False)

# Merge overlapping candidate intervals so the same official region is transferred only once.
intervals=[]
for chrom,d in c.sort_values(['qtl_chr','coloc_region_start']).groupby('qtl_chr'):
    cur=None
    for r in d.itertuples(index=False):
        if cur is None or int(r.coloc_region_start)>cur['end']:
            if cur: intervals.append(cur)
            cur={'chr':chrom,'start':int(r.coloc_region_start),'end':int(r.coloc_region_end),'genes':[r.gene_symbol]}
        else:
            cur['end']=max(cur['end'],int(r.coloc_region_end)); cur['genes'].append(r.gene_symbol)
    if cur: intervals.append(cur)

inventory=[]
for interval_i,iv in enumerate(intervals, start=1):
    region=f"{iv['chr']}:{iv['start']}-{iv['end']}"
    members=c[c.gene_symbol.isin(iv['genes'])]
    gout,gerr=run_tabix(GW,region)
    gd=pd.DataFrame([x.split('\t') for x in gout.splitlines() if x],columns=gcols) if gout.strip() else pd.DataFrame(columns=gcols)
    if len(gd): gd['position']=pd.to_numeric(gd.position,errors='coerce')
    for r in members.itertuples(index=False):
        gf=W/f'{r.gene_symbol}__SCIATICA.tsv.gz'
        gd.to_csv(gf,sep='\t',index=False,compression='gzip')
    # The seven official contexts are independent HTTP range queries.  Query them
    # concurrently within an interval while keeping intervals serial to bound RAM.
    with ThreadPoolExecutor(max_workers=len(contexts)) as pool:
        futures={pool.submit(run_tabix,url,region):(context,dataset,n,url)
                 for context,dataset,n,url in contexts}
        context_results=[]
        for future in as_completed(futures):
            context,dataset,n,url=futures[future]
            out,err=future.result()
            context_results.append((context,dataset,n,url,out,err))
    for context,dataset,n,url,out,err in sorted(context_results, key=lambda x:x[1]):
        qall=pd.DataFrame([x.split('\t') for x in out.splitlines() if x],columns=qcols) if out.strip() else pd.DataFrame(columns=qcols)
        if len(qall): qall['gene_id']=qall.gene_id.astype(str).str.split('.').str[0]
        for r in members.itertuples(index=False):
            gf=W/f'{r.gene_symbol}__SCIATICA.tsv.gz'
            q=qall[qall.gene_id.eq(str(r.ensembl_gene).split('.')[0])].copy() if len(qall) else qall.copy()
            qf=W/f'{r.gene_symbol}__{dataset}.tsv.gz'; q.to_csv(qf,sep='\t',index=False,compression='gzip')
            status='AVAILABLE' if len(q) else ('RETRIEVAL_FAILED' if err else 'NOT_AVAILABLE_FOR_GENE')
            inventory.append(dict(gene_symbol=r.gene_symbol,ensembl_gene=r.ensembl_gene,context=context,dataset_id=dataset,
              qtl_sample_size=n,query_region=region,query_scope='all cis rows returned in official .all file; exact Ensembl gene filter',
              qtl_rows=len(q),gwas_rows=len(gd),qtl_file=str(qf.relative_to(ROOT)),gwas_file=str(gf.relative_to(ROOT)),status=status,error=err))
    print(f'completed interval {interval_i}/{len(intervals)} {region} genes={len(members)}', flush=True)

inv=pd.DataFrame(inventory); inv.to_csv(OUT/'coloc_input_inventory.tsv',sep='\t',index=False)
for f in list(W.glob('*.tsv.gz'))+[OUT/'coloc_input_inventory.tsv']:
    b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print(inv.groupby('status').size().to_dict())
