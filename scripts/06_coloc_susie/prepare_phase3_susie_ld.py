#!/usr/bin/env python3
from pathlib import Path
import subprocess,hashlib
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/06_coloc'; H=ROOT/'work/phase3_final/coloc/harmonized'; L=ROOT/'work/phase3_final/coloc/susie_ld'
L.mkdir(parents=True,exist_ok=True)
r=pd.read_csv(OUT/'coloc_all.tsv',sep='\t'); r=r[r.status.eq('COMPLETED')].sort_values('PP4',ascending=False)
sel=pd.concat([r[r.PP4.ge(.6)],r.groupby('gene_symbol',as_index=False).head(1)],ignore_index=True).drop_duplicates(['gene_symbol','dataset_id'])
plink='plink'; ldroot=ROOT/'work/software/phase3_twas/LDREF'
pairs=[]
for x in sel.itertuples(index=False):
    f=H/f'{x.gene_symbol}__{x.dataset_id}.tsv.gz'
    d=pd.read_csv(f,sep='\t'); chrom=str(d['chr'].dropna().astype(str).iloc[0]).replace('.0','')
    pairs.append(dict(gene_symbol=x.gene_symbol,dataset_id=x.dataset_id,context=x.context,qtl_sample_size=x.qtl_sample_size,PP4_abf=x.PP4,chrom=chrom,harmonized_file=str(f.relative_to(ROOT))))
pm=pd.DataFrame(pairs)
for x in pm.itertuples(index=False):
    d=pd.read_csv(ROOT/x.harmonized_file,sep='\t'); extract=L/f'{x.gene_symbol}__{x.dataset_id}.extract'; extract.write_text('\n'.join(d.rsid.dropna().astype(str).drop_duplicates())+'\n')
    pair=L/f'{x.gene_symbol}__{x.dataset_id}'
    subprocess.run([plink,'--bfile',str(ldroot/f'1000G.EUR.{x.chrom}'),'--extract',str(extract),'--make-bed','--out',str(pair),'--allow-no-sex'],check=True)
    subprocess.run([plink,'--bfile',str(pair),'--r','square','gz','--out',str(pair),'--allow-no-sex'],check=True)
    xdict=x._asdict(); xdict['plink_prefix']=str(pair.relative_to(ROOT)); xdict['n_reference_snps']=sum(1 for _ in open(str(pair)+'.bim')); xdict['ld_source']='official FUSION 1000G EUR per-chromosome reference'
    pairs[pairs.index(next(z for z in pairs if z['gene_symbol']==x.gene_symbol and z['dataset_id']==x.dataset_id))]=xdict
pm=pd.DataFrame(pairs); pm.to_csv(OUT/'susie_ld_manifest.tsv',sep='\t',index=False)
for f in list(L.glob('*'))+[OUT/'susie_ld_manifest.tsv']:
    if f.is_file():
        b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print(pm[['gene_symbol','context','n_reference_snps']].to_string(index=False))
