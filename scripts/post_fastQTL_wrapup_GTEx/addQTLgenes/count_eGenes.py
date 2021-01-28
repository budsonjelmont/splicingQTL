import pandas as pd
import sys

parser = ap.ArgumentParser(description='Read the [pid/sid]-gene_annotations.txt file produced by the bedtools overlaps from count_eGenes.lsf and merge the overlapping gene ID into the QTL results file.')
parser.add_argument('overlapsfile', metavar='overlapsfile', type=str, nargs=1, help='Path to 19-column bedtools intersect output')
parser.add_argument('fqtloutfile', metavar='fqtloutfile', type=str, nargs=1, help='qtls.txt output file to add gene annotation string to')
parser.add_argument('pidorsid', metavar='pidorsid', type=str, nargs=1, help='{\'pid\',\'sid\'} to determine if the overlap is being calculated for SNPs (sid) or phenotypes/introns (pid)')
parser.add_argument('--tag', metavar='tag', type=str, default='', nargs=1, help='(Optional) String to append to the end of the new column name. Default is null string, which will add column "[pid/sid]_ensg" to the qtls.txt file.')
parser.add_argument('--annotate_qtlres', metavar='annotate_qtlres', action='store_true', help='True/False: Should overlapping gene IDs be added to the qtls.txt file? Default: False.')

args = parser.parse_args()

overlapsfile = args.overlapsfile[0]
fqtloutfile = args.fqtloutfile[0]
pidorsid = args.pidorsid[0]
tag = args.tag[0]
annotate_qtlres = args.annotate_qtlres[0]

if pidorsid not in ['pid','sid']:
  print('ERROR: Argument 3 must be one of \'pid\' or \'sid\'')
  exit()
else:
  genecol=pidorsid + '_ensg' + tag

ol = pd.read_csv(overlapsfile, sep='\t', header=None)

# Keep just the 'stable' portion of the gene ID, dropping version info
ol['Gene ID'] = ol[9].str.split('\.',expand=True)[0]

# Read in fastQTL results
fqtlout = pd.read_csv(fqtloutfile,sep='\s+') #'/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.csv'
# Use the phenotype id as the key to grab the (first) overlapping gene ID...
fqtlout.set_index(pidorsid, inplace=True, drop=False, append=False, verify_integrity=False)
fqtlout[genecol] = ol.groupby(3)['Gene ID'].first()

# Print unique genes
print('# of unique ENSGs: ' + str(len(fqtlout[genecol].dropna().unique())))

# Write out the new file with the gene ID included
fqtlout.to_csv(fqtloutfile[:-4] + '+' + genecol + '.txt', sep='\t', index=False, header=True)
