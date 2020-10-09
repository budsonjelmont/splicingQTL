import sys
import pandas as pd
import argparse as ap

parser = ap.ArgumentParser(description='This program is called by find_QTL_overlaps.sh to process the output of the bedtools intersection (-wo) of a BED12 file of IsoSeq Map data with a BED6 file containing SNP data & return the number of overlapping events')
parser.add_argument('overlapbedfile', metavar='overlapbedfile', type=str, nargs=1, help='Path to 18-column bedtools intersect output')
parser.add_argument('backgroundbedfile', metavar='backgroundbedfile', type=str, nargs=1, help='Path BED12 file of all map transcripts used in overlap (-a argument to bedtools)')
parser.add_argument('allqtlsbedfile', metavar='allqtlsbedfile', type=str, nargs=1, help='Path to BED6 file of all QTLs/SNPs (-b argument to bedtools)')
parser.add_argument('statsfile', metavar='statsfile', type=str, nargs=1, help='Output file to concatenate to')
parser.add_argument('comparison', metavar='comparison', type=str, nargs=1, help='Comparison string')

parser.add_argument('-n','--noNMD', dest='noNMD', action='store_true', default=False,help='remove NMD transcripts')

args = parser.parse_args()

overlapbedfile = args.overlapbedfile[0]
backgroundbedfile = args.backgroundbedfile[0]
allqtlsbedfile = args.allqtlsbedfile[0]
statsfile = args.statsfile[0]
comparison = args.comparison[0]

overlapbed = pd.read_csv(overlapbedfile, sep='\t', header=None)
backgroundbed = pd.read_csv(backgroundbedfile, sep='\t', header=None)
allqtlsbed = pd.read_csv(allqtlsbedfile, sep='\t', header=None)
stats = pd.read_csv(statsfile, sep='\t')

stats.set_index('comparison',inplace=True,verify_integrity=True)

def get_map_ids(df,idcol):
  return df[idcol].str.split('|',expand=True)

overlapidsplit = get_map_ids(overlapbed,3) 
backgroundidsplit = get_map_ids(backgroundbed,3)

overlapbed['tscr'] = overlapidsplit[0] 
backgroundbed['tscr'] = backgroundidsplit[0]
overlapbed['gene'] = overlapidsplit[1] 
backgroundbed['gene'] = backgroundidsplit[1]

stats.loc[comparison,'n_overlap_transcripts'] = overlapbed['tscr'].nunique()
stats.loc[comparison,'n_total_transcripts'] = backgroundbed['tscr'].nunique()
stats.loc[comparison,'n_overlap_genes'] = overlapbed['gene'].nunique()
stats.loc[comparison,'n_total_genes'] = backgroundbed['gene'].nunique()
stats.loc[comparison,'n_overlap_snps'] = overlapbed[15].nunique()
stats.loc[comparison,'n_total_snps'] = allqtlsbed[3].nunique()

stats['perc_transcripts'] = stats['n_overlap_transcripts']/stats['n_total_transcripts'] 
stats['perc_genes'] = stats['n_overlap_genes']/stats['n_total_genes'] 
stats['perc_snps'] = stats['n_overlap_snps']/stats['n_total_snps'] 

stats.to_csv(statsfile, sep='\t')
