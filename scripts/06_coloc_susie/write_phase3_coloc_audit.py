#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd

ROOT=Path('.')
OUT=ROOT/'results/phase3_final/06_coloc'
inv=pd.read_csv(OUT/'coloc_input_inventory.tsv',sep='\t')
abf=pd.read_csv(OUT/'coloc_all.tsv',sep='\t')
su=pd.read_csv(OUT/'susie_coloc.tsv',sep='\t')
available=int((inv.status=='AVAILABLE').sum())
missing=int((inv.status=='NOT_AVAILABLE_FOR_GENE').sum())
completed=int((abf.status=='COMPLETED').sum())
strong=int((abf.PP4>=.8).sum())
moderate=int(((abf.PP4>=.6)&(abf.PP4<.8)).sum())
strong_genes=', '.join(sorted(abf.loc[abf.PP4>=.8,'gene_symbol'].unique())) or 'none'
susie_pairs=su[['gene_symbol','dataset_id','status']].drop_duplicates()
susie_signal=int((susie_pairs.status=='COMPLETED').sum())
susie_nosignal=int((susie_pairs.status=='COMPLETED_NO_SIGNAL_PAIR').sum())
susie_strong_genes=', '.join(sorted(su.loc[(su.status=='COMPLETED')&(su['PP.H4.abf']>=.8),'gene_symbol'].unique())) or 'none'
contexts=(inv[['context','dataset_id','qtl_sample_size']].drop_duplicates()
          .sort_values('dataset_id').to_dict('records'))
context_lines='\n'.join(f"- {x['context']}: {x['dataset_id']}, N={int(x['qtl_sample_size'])}" for x in contexts)
text=f'''# QTL data and colocalization audit

## Input contract

- Target universe: 20 non-MHC SCIATICA genetics-only genes selected before spatial or single-cell results were read.
- Scope: official eQTL Catalogue `.all` files queried by tabix over the complete predeclared cis window, then filtered by exact Ensembl gene ID. Significant-only QTL extracts were not used.
- Build: the GWAS and queried eQTL Catalogue positions are GRCh38. Harmonization required shared rsID plus exact REF/effect-allele compatibility; QTL beta was flipped only for a verified allele swap.
- FinnGen case-control input: 28,094 cases and 347,768 controls; no per-variant INFO was invented.
- Coloc priors: p1=1e-4, p2=1e-4, p12=1e-5. PP4 >= 0.80 was called strong and 0.60 <= PP4 < 0.80 moderate.

Predeclared contexts:

{context_lines}

## Availability and results

- Candidate-context contracts: {len(inv)}.
- Full-cis gene/context inputs available: {available}; official interval returned no rows for the exact gene in {missing} contracts.
- ABF coloc completed: {completed}; strong pairs: {strong}; moderate pairs: {moderate}.
- Genes with at least one strong ABF pair: {strong_genes}.
- SuSiE-coloc evaluated {len(susie_pairs)} distinct priority gene/context pairs using signed LD from the official FUSION 1000G EUR reference. Both GWAS and QTL SuSiE fits converged for every reported pair. {susie_signal} pairs produced one or more signal-pair comparisons and {susie_nosignal} produced no signal pair.
- Genes with a SuSiE signal-pair PP4 >= 0.80: {susie_strong_genes}.
- Small negative LD eigenvalues (approximately 1e-5 in magnitude) arise from PLINK text precision; matrices were symmetric and were not silently replaced by simulated LD.

## Interpretation controls

ABF PP4 is not called causality. Multi-signal SuSiE-coloc is the robustness layer: ABF-strong results without a SuSiE signal pair remain provisional rather than being presented as multi-signal-confirmed. The chromosome 3 and chromosome 12 target clusters contain several nearby genes, so competing-gene evidence is shown together and no lead gene is selected by PP4 alone. Lead-SNP direction is reported as allele-aligned concordant or discordant; it is not interpreted as whether increased expression is biologically beneficial without a validated expression model.
'''
f=OUT/'qtl_data_audit.md'; f.write_text(text)
for p in [OUT/'coloc_all.tsv',OUT/'coloc_strong.tsv',OUT/'coloc_harmonization.tsv',OUT/'susie_coloc.tsv',OUT/'susie_ld_manifest.tsv',OUT/'coloc_input_inventory.tsv',f]:
    b=p.read_bytes()
    Path(str(p)+'.md5').write_text(hashlib.md5(b).hexdigest()+'  '+p.name+'\n')
    Path(str(p)+'.sha256').write_text(hashlib.sha256(b).hexdigest()+'  '+p.name+'\n')
print({'available':available,'missing':missing,'ABF_completed':completed,'ABF_strong':strong,'ABF_moderate':moderate,'SuSiE_distinct_pairs':len(susie_pairs),'SuSiE_signal_pairs':susie_signal})
