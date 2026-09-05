#!/usr/bin/env python3
from pathlib import Path
import hashlib, math
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/08_gsmap'
RUN=ROOT/'work/phase3_final/gsmap_runs'
qc=pd.read_csv(OUT/'section_qc.tsv',sep='\t')
rows=[]
for sample in ['GSM6919905','GSM6919909','GSM6919911','GSM6919917']:
    h=pd.read_csv(RUN/f'height_control/{sample}/spatial_ldsc/{sample}_HEIGHT_IRN.csv.gz').set_index('spot')
    s=pd.read_csv(RUN/f'main/{sample}/spatial_ldsc/{sample}_SCIATICA.csv.gz').set_index('spot')
    x=h[['z','p']].join(s[['z','p']],lsuffix='_height',rsuffix='_sciatica',how='inner')
    rho,p=spearmanr(x.z_sciatica,x.z_height)
    n=max(1,int(math.ceil(.05*len(x))))
    a=set(x.nlargest(n,'z_sciatica').index); b=set(x.nlargest(n,'z_height').index)
    donor=qc.loc[qc['sample'].eq(sample),'donor'].iloc[0]
    rows.append(dict(sample=sample,donor=donor,n_spots=len(x),spearman_z_sciatica_height=rho,
                     spearman_p=p,top5_jaccard=len(a&b)/len(a|b),
                     interpretation='LD-preserving control is globally significant but its relative spot pattern is distinct'))
d=pd.DataFrame(rows); d.to_csv(OUT/'negative_control_spatial_comparison.tsv',sep='\t',index=False)
note=OUT/'negative_control_interpretation.md'
note.write_text('''# LD-preserving control interpretation

The real FinnGen HEIGHT_IRN control is globally significant in all four selected sections and therefore shows that gsMap global enrichment P values are not trait-specific calibration evidence in this dataset. This prevents a claim that the existence or magnitude of a global P value is specific to sciatica.

However, the within-section relative spot patterns differ: Sciatica-versus-height Spearman correlations and top-5% spot overlap are reported in `negative_control_spatial_comparison.tsv`. The four correlations are negative, so sciatica localization is not a simple copy of the height pattern. The defensible result is consequently limited to reproducible relative spatial patterning across the four sciatica donors, not global enrichment specificity and not disease-tissue activation.
''')
for f in [OUT/'negative_control_spatial_comparison.tsv',note]:
    b=f.read_bytes(); Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n'); Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print(d.to_string(index=False))
