#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd
ROOT=Path('.')
O=ROOT/'results/phase3_final/08_gsmap'; rows=[]
for trait in ['SCIATICA','LDH']:
 a=pd.read_csv(ROOT/f'work/phase3_final/gsmap_runs/main/GSM6919905/spatial_ldsc/GSM6919905_{trait}.csv.gz')
 b=pd.read_csv(ROOT/f'work/phase3_final/gsmap_runs/no_mhc/GSM6919905/spatial_ldsc/GSM6919905_{trait}_noMHC.csv.gz')
 d=a.merge(b,on='spot',suffixes=('_main','_noMHC'))
 rows.append(dict(sample='GSM6919905',trait=trait,n_spots=len(d),pearson_z=d.z_main.corr(d.z_noMHC),max_abs_z_delta=(d.z_main-d.z_noMHC).abs().max(),max_abs_p_delta=(d.p_main-d.p_noMHC).abs().max(),interpretation='identical because distributed gsMap LD weights are HM3 no-HLA'))
f=O/'mhc_sensitivity_GSM6919905.tsv';pd.DataFrame(rows).to_csv(f,sep='\t',index=False)
n=O/'mhc_sensitivity_validation_note.md';n.write_text('''# gsMap MHC sensitivity validation\n\nThe first direct rerun command exited 0 but found zero retained chunk files and produced no sensitivity CSV; product validation therefore marked it unusable. A complete independent quick-mode rerun was then performed on GSM6919905 with explicit MHC-removed sumstats. Both SCIATICA and LDH spot statistics were exactly identical to the main run. This is expected, not an independent numerical coincidence: the distributed gsMap regression-weight reference is `weights_hm3_no_hla`, so MHC variants never enter the final SNP intersection in either input.\n''')
for x in [f,n]:
 b=x.read_bytes();Path(str(x)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+x.name+'\n');Path(str(x)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+x.name+'\n')
print(pd.DataFrame(rows).to_string(index=False))
