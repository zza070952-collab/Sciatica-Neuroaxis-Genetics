#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd

ROOT = Path('.')
OUT = ROOT / 'results/phase3_final/05_twas'
d = pd.read_csv(OUT / 'twas_all_results.tsv', sep='\t')
j = pd.read_csv(OUT / 'conditional_joint_results.tsv', sep='\t')

rows = []
for trait, x in d.groupby('trait'):
    hit = x[x.fdr < .05]
    multi = hit[hit.n_snps_used >= 2]
    bonf = x[x.pvalue < x.bonferroni_trait]
    jx = j[(j.trait == trait) & (j.joint_status == 'INCLUDED')] if 'trait' in j else j.iloc[0:0]
    rows.append({
        'trait': trait, 'tested_models': len(x), 'tested_unique_genes': x.ensembl_gene.nunique(),
        'traitwide_FDR_models': len(hit), 'traitwide_FDR_unique_genes': hit.ensembl_gene.nunique(),
        'multi_SNP_traitwide_FDR_models': len(multi), 'multi_SNP_traitwide_FDR_unique_genes': multi.ensembl_gene.nunique(),
        'Bonferroni_models': len(bonf), 'Bonferroni_unique_genes': bonf.ensembl_gene.nunique(),
        'joint_included_models': len(jx), 'joint_included_unique_genes': jx.ID.astype(str).str.split('.').str[0].nunique() if len(jx) else 0
    })
s = pd.DataFrame(rows)
sci = set(d[(d.trait == 'SCIATICA') & (d.fdr < .05) & (d.n_snps_used >= 2)].ensembl_gene)
ldh = set(d[(d.trait == 'LDH') & (d.fdr < .05) & (d.n_snps_used >= 2)].ensembl_gene)
s['cross_trait_multi_SNP_FDR_gene_overlap'] = len(sci & ldh)
s.to_csv(OUT / 'twas_summary.tsv', sep='\t', index=False)

top = d[(d.trait == 'SCIATICA') & (d.fdr < .05) & (d.n_snps_used >= 2)].sort_values(['pvalue', 'fdr'])
top[['gene_name','ensembl_gene','tissue','zscore','pvalue','fdr','n_snps_used','MODEL']].head(30).to_csv(
    OUT / 'twas_top_multisnp_sciatica.tsv', sep='\t', index=False)

report = f'''# Formal FUSION TWAS summary

The five predeclared GTEx v8 EUR tissues were analyzed with FUSION and matched 1000 Genomes EUR LD. Each trait tested {int(s.loc[s.trait.eq('SCIATICA'),'tested_models'].iloc[0]):,} tissue-gene models. Across the five tissues, sciatica had {int(s.loc[s.trait.eq('SCIATICA'),'traitwide_FDR_models'].iloc[0]):,} trait-wide FDR-significant models ({int(s.loc[s.trait.eq('SCIATICA'),'multi_SNP_traitwide_FDR_models'].iloc[0]):,} multi-SNP), while LDH had {int(s.loc[s.trait.eq('LDH'),'traitwide_FDR_models'].iloc[0]):,} ({int(s.loc[s.trait.eq('LDH'),'multi_SNP_traitwide_FDR_models'].iloc[0]):,} multi-SNP). The multi-SNP FDR gene sets overlap at {len(sci & ldh)} stable Ensembl genes; this is cross-trait concordance, not independence.

FUSION post-processing used a trait-wide Bonferroni Z threshold and the same EUR LD. It retained {int(s.loc[s.trait.eq('SCIATICA'),'joint_included_models'].iloc[0])} sciatica and {int(s.loc[s.trait.eq('LDH'),'joint_included_models'].iloc[0])} LDH tissue-gene models. One-SNP models remain visible but cannot independently qualify a candidate. Formal targeted full-cis colocalization, not TWAS significance alone, defines stronger regulatory support.
'''
(OUT / 'twas_summary.md').write_text(report)

for path in [OUT/'twas_summary.tsv', OUT/'twas_top_multisnp_sciatica.tsv', OUT/'twas_summary.md']:
    data = path.read_bytes()
    Path(str(path)+'.md5').write_text(hashlib.md5(data).hexdigest()+'  '+path.name+'\n')
    Path(str(path)+'.sha256').write_text(hashlib.sha256(data).hexdigest()+'  '+path.name+'\n')
print(s.to_string(index=False))
print(top[['gene_name','tissue','zscore','pvalue','fdr','n_snps_used','MODEL']].head(20).to_string(index=False))
