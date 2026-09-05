#!/usr/bin/env python3
from pathlib import Path
from statistics import NormalDist
import hashlib
import re
import numpy as np
import pandas as pd

ROOT = Path('.')
OUT = ROOT / 'results/phase3_final/05_twas'
RAW = OUT / 'fusion_raw'
JOINT = OUT / 'fusion_joint'
GTF = ROOT / 'work/phase3_final/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'
TISSUES = ['Nerve_Tibial', 'Brain_Spinal_cord_cervical_c-1', 'Whole_Blood', 'Muscle_Skeletal', 'Cells_Cultured_fibroblasts']
TRAITS = ['SCIATICA', 'LDH']

def bh(values):
    p = np.asarray(values, dtype=float)
    n = len(p)
    order = np.argsort(np.where(np.isfinite(p), p, 1.0))
    q = np.ones(n, dtype=float)
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, 1))):
        val = p[idx] if np.isfinite(p[idx]) else 1.0
        running = min(running, val * n / rank)
        q[idx] = running
    return q

def checksum(path):
    data = path.read_bytes()
    Path(str(path) + '.md5').write_text(hashlib.md5(data).hexdigest() + '  ' + path.name + '\n')
    Path(str(path) + '.sha256').write_text(hashlib.sha256(data).hexdigest() + '  ' + path.name + '\n')

gene_map = {}
with GTF.open(errors='replace') as handle:
    for line in handle:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) < 9 or fields[2] != 'gene':
            continue
        attrs = dict(re.findall(r'(\S+) "([^"]+)"', fields[8]))
        if 'gene_id' in attrs and 'gene_name' in attrs:
            gene_map[attrs['gene_id'].split('.')[0]] = attrs['gene_name']

rows = []
product_rows = []
combined_dir = OUT / 'fusion_combined_input'
combined_dir.mkdir(parents=True, exist_ok=True)
for trait in TRAITS:
    combined_parts = []
    for tissue in TISSUES:
        for chrom in range(1, 23):
            path = RAW / trait / tissue / f'chr{chrom}.dat'
            if not path.exists() or path.stat().st_size == 0:
                raise RuntimeError(f'missing FUSION product: {path}')
            frame = pd.read_csv(path, sep='\t')
            required = {'ID', 'CHR', 'P0', 'P1', 'NSNP', 'NWGT', 'MODEL', 'MODELCV.R2', 'MODELCV.PV', 'TWAS.Z', 'TWAS.P'}
            if not required.issubset(frame.columns):
                raise RuntimeError(f'invalid FUSION columns: {path}')
            product_rows.append({'trait': trait, 'tissue': tissue, 'chromosome': chrom, 'rows': len(frame), 'status': 'PASS'})
            frame.insert(0, 'trait', trait)
            frame.insert(1, 'tissue', tissue)
            rows.append(frame)
            combined_parts.append(frame.drop(columns=['trait', 'tissue']))
    combined = pd.concat(combined_parts, ignore_index=True)
    # FUSION.post_process.R uses read.table with whitespace separation; explicit
    # NA tokens are required because empty tab fields would be collapsed.
    combined.to_csv(combined_dir / f'{trait}.dat', sep='\t', index=False, na_rep='NA')

d = pd.concat(rows, ignore_index=True)
d['gene'] = d['ID'].astype(str)
d['ensembl_gene'] = d['gene'].str.split('.').str[0]
d['gene_name'] = d['ensembl_gene'].map(gene_map).fillna(d['ensembl_gene'])
d['pvalue'] = pd.to_numeric(d['TWAS.P'], errors='coerce')
d['zscore'] = pd.to_numeric(d['TWAS.Z'], errors='coerce')
d['n_snps_used'] = pd.to_numeric(d['NWGT'], errors='coerce').astype('Int64')
d['n_locus_snps'] = pd.to_numeric(d['NSNP'], errors='coerce').astype('Int64')
d['fdr_tissue'] = d.groupby(['trait', 'tissue'])['pvalue'].transform(bh)
d['fdr'] = d.groupby('trait')['pvalue'].transform(bh)
d['one_snp_model'] = d['n_snps_used'].eq(1)
d['method'] = 'FUSION'
d['weights'] = 'GTEx_v8_EUR'
d['ld_reference'] = '1000G_EUR'
d['calibration_status'] = 'VALIDATED_FUSION_Z_WITH_MATCHED_EUR_LD'
d['bonferroni_trait'] = d.groupby('trait')['pvalue'].transform(lambda x: 0.05 / len(x))

cols = ['trait', 'tissue', 'gene', 'ensembl_gene', 'gene_name', 'CHR', 'P0', 'P1', 'zscore', 'pvalue',
        'fdr', 'fdr_tissue', 'bonferroni_trait', 'n_snps_used', 'n_locus_snps', 'MODEL', 'MODELCV.R2',
        'MODELCV.PV', 'HSQ', 'BEST.GWAS.ID', 'BEST.GWAS.Z', 'EQTL.ID', 'EQTL.R2', 'EQTL.Z',
        'EQTL.GWAS.Z', 'one_snp_model', 'method', 'weights', 'ld_reference', 'calibration_status', 'FILE', 'PANEL']
d[cols].to_csv(OUT / 'twas_all_results.tsv', sep='\t', index=False)
d.loc[d['fdr'] < 0.05, cols].to_csv(OUT / 'twas_FDR_results.tsv', sep='\t', index=False)

manifest = []
for (trait, tissue), frame in d.groupby(['trait', 'tissue']):
    manifest.append({
        'trait': trait, 'tissue': tissue, 'method': 'FUSION', 'weights': 'GTEx_v8_EUR',
        'tested_models': len(frame), 'tested_genes': frame['ensembl_gene'].nunique(),
        'trait_wide_FDR_lt_0.05': int((frame['fdr'] < 0.05).sum()),
        'tissue_wide_FDR_lt_0.05': int((frame['fdr_tissue'] < 0.05).sum()),
        'multi_SNP_trait_wide_FDR_lt_0.05': int(((frame['fdr'] < 0.05) & (frame['n_snps_used'] >= 2)).sum()),
        'build': 'GRCh37', 'ancestry': 'EUR', 'status': 'COMPLETE'
    })
pd.DataFrame(manifest).to_csv(OUT / 'tissue_manifest.tsv', sep='\t', index=False)
pd.DataFrame(product_rows).to_csv(OUT / 'fusion_product_validation.tsv', sep='\t', index=False)

comparison = d.pivot_table(index=['gene', 'ensembl_gene', 'gene_name'], columns=['trait', 'tissue'], values='pvalue').reset_index()
comparison.columns = ['|'.join(str(y) for y in x if str(y)) if isinstance(x, tuple) else x for x in comparison.columns]
comparison.to_csv(OUT / 'tissue_comparison.tsv', sep='\t', index=False)

thresholds = []
for trait, frame in d.groupby('trait'):
    alpha = 0.05 / len(frame)
    thresholds.append({'trait': trait, 'tested_models': len(frame), 'bonferroni_p': alpha,
                       'two_sided_abs_z': NormalDist().inv_cdf(1 - alpha / 2)})
pd.DataFrame(thresholds).to_csv(OUT / 'fusion_bonferroni_thresholds.tsv', sep='\t', index=False)

joint_rows = []
if JOINT.exists():
    for path in sorted(JOINT.glob('*/*.joint_included.dat')):
        trait = path.parent.name
        match = re.search(r'chr(\d+)', path.name)
        chrom = int(match.group(1)) if match else np.nan
        try:
            frame = pd.read_csv(path, sep='\t')
        except pd.errors.EmptyDataError:
            continue
        if len(frame):
            frame.insert(0, 'trait', trait)
            frame.insert(1, 'chromosome', chrom)
            frame.insert(2, 'joint_status', 'INCLUDED')
            joint_rows.append(frame)
    for path in sorted(JOINT.glob('*/*.joint_dropped.dat')):
        trait = path.parent.name
        match = re.search(r'chr(\d+)', path.name)
        chrom = int(match.group(1)) if match else np.nan
        try:
            frame = pd.read_csv(path, sep='\t')
        except pd.errors.EmptyDataError:
            continue
        if len(frame):
            frame.insert(0, 'trait', trait)
            frame.insert(1, 'chromosome', chrom)
            frame.insert(2, 'joint_status', 'DROPPED')
            joint_rows.append(frame)
if joint_rows:
    pd.concat(joint_rows, ignore_index=True).to_csv(OUT / 'conditional_joint_results.tsv', sep='\t', index=False)
else:
    pd.DataFrame([{'analysis': 'conditional_joint', 'status': 'PENDING_FUSION_POST_PROCESS'}]).to_csv(OUT / 'conditional_joint_results.tsv', sep='\t', index=False)

audit = '''# Phase 3 TWAS method audit

The formal TWAS uses official FUSION GTEx v8 EUR weights in exactly five predeclared tissues and the 1000 Genomes EUR LD reference. FinnGen LDSC-format HapMap3 inputs were restricted to complete SNP/A1/A2/Z/N rows; incomplete rows were removed and counted. Trait-wide Benjamini-Hochberg FDR (`fdr`) is calculated across every model in all five predeclared tissues for each trait, preventing tissue cherry-picking. Tissue-specific FDR is retained as a secondary field. One-SNP models are flagged and cannot independently qualify a prioritized candidate.


FUSION conditional/joint analysis uses the same five-tissue models and matched EUR LD. TWAS is regulatory association evidence, not proof that a gene is causal; targeted full-cis colocalization remains required for stronger regulatory support.
'''
(OUT / 'twas_methods_audit.md').write_text(audit)

for path in [OUT / 'twas_all_results.tsv', OUT / 'twas_FDR_results.tsv', OUT / 'tissue_manifest.tsv',
             OUT / 'fusion_product_validation.tsv', OUT / 'tissue_comparison.tsv', OUT / 'conditional_joint_results.tsv',
             OUT / 'fusion_bonferroni_thresholds.tsv', OUT / 'twas_methods_audit.md']:
    checksum(path)
for path in combined_dir.glob('*.dat'):
    checksum(path)

print({'models': len(d), 'trait_wide_FDR_hits': int((d['fdr'] < 0.05).sum()),
       'multi_SNP_hits': int(((d['fdr'] < 0.05) & (d['n_snps_used'] >= 2)).sum()),
       'joint_rows': sum(len(x) for x in joint_rows) if joint_rows else 0})
