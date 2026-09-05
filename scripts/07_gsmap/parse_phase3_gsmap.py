#!/usr/bin/env python3
from pathlib import Path
import hashlib
import math
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import norm, chi2, spearmanr

ROOT = Path('.')
OUT = ROOT / 'results/phase3_final/08_gsmap'
RUN = ROOT / 'work/phase3_final/gsmap_runs/main'
H5 = ROOT / 'work/phase3_final/gsmap_input_pca'
qc = pd.read_csv(OUT / 'section_qc.tsv', sep='\t')

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p); q = np.empty(len(p)); prev = 1.0
    for rank, idx in reversed(list(enumerate(o, 1))):
        prev = min(prev, p[idx] * len(p) / rank); q[idx] = prev
    return q

def acat(p):
    p = np.clip(np.asarray(p, float), 1e-15, 1-1e-15)
    t = np.mean(np.tan((0.5-p)*np.pi))
    return float(np.clip(0.5-np.arctan(t)/np.pi, 0, 1))

spot_parts=[]; sections=[]
for r in qc.itertuples(index=False):
    a = ad.read_h5ad(H5 / f'{r.sample}.h5ad', backed='r')
    xy = pd.DataFrame(np.asarray(a.obsm['spatial']), index=a.obs_names,
                      columns=['spatial_x','spatial_y'])
    for trait in ['SCIATICA','LDH']:
        f = RUN / r.sample / 'spatial_ldsc' / f'{r.sample}_{trait}.csv.gz'
        if not f.exists():
            raise FileNotFoundError(f)
        d = pd.read_csv(f).set_index('spot').join(xy, how='left').reset_index()
        if d[['spatial_x','spatial_y']].isna().any().any():
            raise RuntimeError(f'coordinate mismatch: {r.sample} {trait}')
        d['q_bh_section'] = bh(d.p)
        d.insert(0,'trait',trait); d.insert(0,'donor',r.donor); d.insert(0,'section',r.section); d.insert(0,'sample',r.sample)
        spot_parts.append(d)
        cf = RUN / r.sample / 'cauchy_combination' / f'{r.sample}_{trait}.Cauchy.csv.gz'
        c = pd.read_csv(cf).iloc[0]
        cutoff = np.nanquantile(d.z, .95)
        sections.append(dict(sample=r.sample,donor=r.donor,section=r.section,trait=trait,
          n_spots=len(d),global_p_cauchy=float(c.p_cauchy),global_p_median=float(c.p_median),
          median_z=float(d.z.median()),top5_median_z=float(d.loc[d.z>=cutoff,'z'].median()),
          max_z=float(d.z.max()),min_p=float(d.p.min()),n_fdr_0_05=int((d.q_bh_section<.05).sum()),
          fraction_p_lt_0_05=float((d.p<.05).mean())))

spots=pd.concat(spot_parts,ignore_index=True)
sec=pd.DataFrame(sections)
spots.to_csv(OUT/'spot_level_results.tsv.gz',sep='\t',index=False,compression='gzip')
sec.to_csv(OUT/'section_level_results.tsv',sep='\t',index=False)

donors=[]
for (donor,trait),d in sec.groupby(['donor','trait']):
    p=acat(d.global_p_cauchy)
    donors.append(dict(donor=donor,trait=trait,n_sections=len(d),n_spots=int(d.n_spots.sum()),
      p_acat_sections=p,z_from_one_sided_p=float(norm.isf(max(p,1e-300))),
      median_section_z=float(d.median_z.median()),median_top5_z=float(d.top5_median_z.median()),
      all_sections_nominal=bool((d.global_p_cauchy<.05).all()),
      n_sections_global_p_lt_0_05=int((d.global_p_cauchy<.05).sum())))
don=pd.DataFrame(donors)
don.to_csv(OUT/'donor_level_results.tsv',sep='\t',index=False)

meta=[]; loo=[]
for trait,d in don.groupby('trait'):
    p_acat=acat(d.p_acat_sections)
    fish=-2*np.log(np.clip(d.p_acat_sections,1e-300,1)).sum()
    p_fisher=float(chi2.sf(fish,2*len(d)))
    meta.append(dict(trait=trait,n_donors=len(d),p_acat_donors=p_acat,p_fisher_donors=p_fisher,
      z_acat=float(norm.isf(max(p_acat,1e-300))),n_donors_p_lt_0_05=int((d.p_acat_sections<.05).sum()),
      all_donors_nominal=bool((d.p_acat_sections<.05).all()),
      donor_z_median=float(d.z_from_one_sided_p.median()),donor_z_sd=float(d.z_from_one_sided_p.std(ddof=1))))
    for leave in d.donor:
        x=d[d.donor.ne(leave)]
        p=acat(x.p_acat_sections)
        loo.append(dict(trait=trait,left_out_donor=leave,n_donors=len(x),p_acat=p,
                        z_acat=float(norm.isf(max(p,1e-300))),all_remaining_nominal=bool((x.p_acat_sections<.05).all())))
pd.DataFrame(meta).to_csv(OUT/'donor_meta_results.tsv',sep='\t',index=False)
pd.DataFrame(loo).to_csv(OUT/'leave_one_donor_out.tsv',sep='\t',index=False)

comparisons=[]
for sample,d in spots.groupby('sample'):
    a=d[d.trait.eq('SCIATICA')].set_index('spot'); b=d[d.trait.eq('LDH')].set_index('spot')
    common=a.index.intersection(b.index); rho,pv=spearmanr(a.loc[common,'z'],b.loc[common,'z'])
    na=max(1,int(math.ceil(.05*len(common))))
    sa=set(a.loc[common].nlargest(na,'z').index); sb=set(b.loc[common].nlargest(na,'z').index)
    comparisons.append(dict(level='section',unit=sample,donor=d.donor.iloc[0],n_spots=len(common),
      spearman_z=float(rho),spearman_p=float(pv),top5_jaccard=len(sa&sb)/len(sa|sb),
      interpretation='spatial concordance; not disease differential expression'))
for donor,d in pd.DataFrame(comparisons).groupby('donor'):
    comparisons.append(dict(level='donor_summary',unit=donor,donor=donor,n_spots=int(d.n_spots.sum()),
      spearman_z=float(d.spearman_z.median()),spearman_p=np.nan,top5_jaccard=float(d.top5_jaccard.median()),
      interpretation='median across sections; sciatica remains the nomination phenotype'))
pd.DataFrame(comparisons).to_csv(OUT/'sciatica_LDH_spatial_comparison.tsv',sep='\t',index=False)

# Real-trait LD-preserving control, one QC-selected section per donor.
neg=[]
for sample in ['GSM6919905','GSM6919909','GSM6919911','GSM6919917']:
    f=ROOT/f'work/phase3_final/gsmap_runs/height_control/{sample}/cauchy_combination/{sample}_HEIGHT_IRN.Cauchy.csv.gz'
    if f.exists():
        x=pd.read_csv(f).iloc[0]; donor=qc.loc[qc['sample'].eq(sample),'donor'].iloc[0]
        neg.append(dict(sample=sample,donor=donor,control_trait='HEIGHT_IRN',global_p_cauchy=float(x.p_cauchy),global_p_median=float(x.p_median),control_role='LD-preserving unrelated quantitative-trait control'))
if neg:
    nd=pd.DataFrame(neg); nd.loc[len(nd)]={'sample':'DONOR_META','donor':'4_donors','control_trait':'HEIGHT_IRN','global_p_cauchy':acat(nd.global_p_cauchy),'global_p_median':float(nd.global_p_median.median()),'control_role':'ACAT across one QC-selected section per donor'}
else:
    nd=pd.DataFrame([dict(sample='NOT_COMPLETED',donor='NA',control_trait='HEIGHT_IRN',global_p_cauchy=np.nan,global_p_median=np.nan,control_role='control products unavailable')])
nd.to_csv(OUT/'negative_control_results.tsv',sep='\t',index=False)

audit='''# gsMap methods audit

- Dataset: GSE222322, 20 sections from four adult donors; each section was processed separately and donor was the biological replicate.
- Coordinates: official barcode-matched coordinates that passed the earlier >=99% gate; no sections were stitched into a pseudo-spatial object.
- Software: gsMap 1.73.7. A verified read-only algorithm-specific environment was used; no previous results or candidates were imported.
- Stable execution mode: the official `latent_representation` option with independently generated PCA50 (HVG=3000, library-size 10,000, log1p, seed 1907). The initial GNN run was stopped after excessive runtime and is retained as a failed attempt, not a result.
- GWAS input: FinnGen R13 common biallelic high-quality release SNP set. The distributed gsMap reference is HapMap3/HLA-excluded, so the primary intersection is already MHC-excluded at the LD-reference stage.
- Spot P values are one-sided positive enrichment tests. BH values are calculated within section and trait. Top-5% spots are descriptive localization neighborhoods when section-wise FDR is not reached.
- Within-donor section P values and across-donor P values were combined using the Cauchy combination test, which is robust to dependence. Fisher combination is reported as a secondary donor-level summary.
- The public matrices contain normal adult tissue. These analyses cannot establish case-control expression, disease-specific activation, or treatment response.
- The LD-preserving control is the real FinnGen R13 inverse-normalized adult-height GWAS, run on one QC-selected section per donor. It preserves a real polygenic/LD architecture and is not a permuted null; it is therefore a calibration control, not proof of biological specificity.
- Official spot-level dorsal/ventral or horn labels were not present in the processed object. Named anatomical subregions are therefore not invented; interpretation is limited to observed spatial neighborhoods and reference-cell localization.
- Donor19 matrices were generated with a reference explicitly labelled `GRCh38_No_Mt_No_LncRNA_No_Pseudogene`; mitochondrial percentage is structurally unavailable (reported as 0), not evidence of zero mitochondrial transcription.
'''
(OUT/'gsmap_methods_audit.md').write_text(audit)

for f in [OUT/'spot_level_results.tsv.gz',OUT/'section_level_results.tsv',OUT/'donor_level_results.tsv',OUT/'donor_meta_results.tsv',OUT/'leave_one_donor_out.tsv',OUT/'sciatica_LDH_spatial_comparison.tsv',OUT/'negative_control_results.tsv',OUT/'gsmap_methods_audit.md']:
    b=f.read_bytes()
    Path(str(f)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+f.name+'\n')
    Path(str(f)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+f.name+'\n')
print({'spots':len(spots),'sections':len(sec),'donors':len(don)})
