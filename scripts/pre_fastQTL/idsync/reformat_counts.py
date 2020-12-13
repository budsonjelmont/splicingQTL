import argparse as ap
import pandas as pd

# Argument parser setup
parser = ap.ArgumentParser(description='Process leafcutter output file, discarding samples that are not included in the \'PSICountID\' column of the genotype-phenotype synchronized metadata file')
parser.add_argument('incountsfile', metavar='incountsfile', type=str, nargs=1, help='Path to the leafcutter output file (*_perind_numers.counts or *_perind.counts)')
parser.add_argument('metafile', metavar='metafile', type=str, nargs=1, help='Path to the synchronized metadata file containing only the samples to be used in the QTL detection. Sample IDs are read from the column \'PSICountID\'')
parser.add_argument('outfile', metavar='outfile', type=str, nargs=1, help='Path to write output file')

args = parser.parse_args()
#args = parser.parse_args('/sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts /sc/arion/projects/EPIASD/splicingQTL/data/metadata/meta_matchedIDs.csv /sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync_test'.split(' '))

incountsfile=args.incountsfile[0]
metafile=args.metafile[0]
outfile=args.outfile[0]

incounts = pd.read_csv(incountsfile,sep=' ')
meta = pd.read_csv(metafile)

# Make list of columns to keep
keepme = ['chrom']
keepme.extend(meta['PSIcountID'].values)

# Write only included rows
incounts.to_csv(outfile, sep=' ', header=True, index=False, columns=keepme)
