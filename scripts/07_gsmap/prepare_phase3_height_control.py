#!/usr/bin/env python3
from pathlib import Path
import gzip,csv,hashlib,math
ROOT=Path('.')
src=Path('data/gwas/finngen_R13_HEIGHT_IRN.gz')
out=ROOT/'work/phase3_final/gwas_formats/HEIGHT_IRN.gsmap.sumstats.gz'
out.parent.mkdir(parents=True,exist_ok=True); counts={'raw':0,'kept':0,'non_snv':0,'maf':0,'ambiguous':0,'rsid':0}
with gzip.open(src,'rt') as fi,gzip.open(str(out)+'.tmp','wt',newline='') as fo:
    r=csv.DictReader(fi,delimiter='\t'); w=csv.writer(fo,delimiter='\t',lineterminator='\n'); w.writerow(['SNP','A1','A2','Z','N'])
    for x in r:
        counts['raw']+=1; ref=x['ref'].upper(); alt=x['alt'].upper(); rs=x['rsids']; af=float(x['af_alt'])
        if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT': counts['non_snv']+=1; continue
        if not (0.01<=af<=0.99): counts['maf']+=1; continue
        if {ref,alt} in [set('AT'),set('CG')] and .4<=af<=.6: counts['ambiguous']+=1; continue
        if not rs.startswith('rs') or ',' in rs or ';' in rs: counts['rsid']+=1; continue
        se=float(x['sebeta']); beta=float(x['beta']);
        if not math.isfinite(se) or se<=0 or not math.isfinite(beta): continue
        w.writerow([rs,alt,ref,beta/se,364515]); counts['kept']+=1
Path(str(out)+'.tmp').replace(out)
(ROOT/'results/phase3_final/08_gsmap/height_control_qc.tsv').write_text('metric\tvalue\n'+''.join(f'{k}\t{v}\n' for k,v in counts.items())+'control_role\tLD-preserving unrelated quantitative-trait control; not a permuted null\n')
for f in [src,out,ROOT/'results/phase3_final/08_gsmap/height_control_qc.tsv']:
    m=hashlib.md5(); s=hashlib.sha256()
    with f.open('rb') as fh:
        for b in iter(lambda:fh.read(8*1024*1024),b''): m.update(b);s.update(b)
    Path(str(f)+'.md5').write_text(m.hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(s.hexdigest()+'  '+f.name+'\n')
print(counts)
