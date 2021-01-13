#!/usr/bin/env python3
# Adapted from https://github.com/francois-a/fastqtl/blob/master/python/annotate_outputs.py
# Calls QTL based on nominal p-value, and npval threshold calculated from permutations
# Process by chunks

import argparse
import numpy as np
import pandas as pd
import os
import gzip
from datetime import datetime
import subprocess
import io

parser = argparse.ArgumentParser(description='Filter significant SNP-gene pairs from FastQTL results using FDR cutoff')
parser.add_argument('--egenes_file', help='eGenes file') # Permutation pass results
parser.add_argument('--nominal_results', help='FastQTL output from nominal pass')
parser.add_argument('-o', '--output_dir', help='Output directory')
args = parser.parse_args()
#args = parser.parse_args('--egenes_file /sc/arion/projects/EPIASD/splicingQTL/output/sqtls/minCovars+seqPC9/4genoPCs/35HCPs/permute/wrapup/chrAll_combined.FDR05 --nominal_results /sc/arion/projects/EPIASD/splicingQTL/output/sqtls/minCovars+seqPC9/4genoPCs/35HCPs/nominal/wrapup/chrAll_combined --output_dir /sc/arion/projects/EPIASD/splicingQTL/output/sqtls/minCovars+seqPC9/4genoPCs/35HCPs/qtls'.split(' '))

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Filtering significant variant-gene pairs', flush=True)
# eGenes
egene_df = pd.read_csv(args.egenes_file, sep='\t', index_col=0)[['pval_nominal_threshold', 'npval', 'bpval']].copy()
egene_df.rename(columns={'npval': 'min_pval_nominal'}, inplace=True)
egene_ids = set(egene_df.index)
threshold_dict = egene_df['pval_nominal_threshold'].to_dict()
# process by chunks to reduce memory usage
signif_df = []
mask = []
for i,chunk in enumerate(pd.read_csv(args.nominal_results, sep='\s+', iterator=True, chunksize=1000000, index_col=1, header=None, names=['pid','sid','dist','npval','slope'])):
    chunk['sid'] = chunk.index # Added so that sid would be retained in output
    chunk = chunk[chunk['pid'].isin(egene_ids)]
    m = chunk['npval'].astype(float)<chunk['pid'].apply(lambda x: threshold_dict[x])
    signif_df.append(chunk[m])
    mask.append(m)
    print('Chunks processed: {0:d}'.format(i+1), end='\r', flush=True)
signif_df = pd.concat(signif_df, axis=0)
signif_df = signif_df.merge(egene_df, left_on='pid', right_index=True)
outname = os.path.join(args.output_dir, 'qtl_chunk_concat.txt')
with open(outname, 'wt') as f:
    signif_df.to_csv(f, sep='\t', float_format='%.6g',index=False) # Remove index in output since we have it already in the sid column

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Completed', flush=True)
