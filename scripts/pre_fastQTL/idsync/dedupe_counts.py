import argparse as ap
import pandas as pd

# Argument parser setup
parser = ap.ArgumentParser(description='Process leafcutter output file, discarding samples that are included in the list passed in the 2nd argument)')
parser.add_argument('incountsfile', metavar='incountsfile', type=str, nargs=1, help='Path to the leafcutter output file (*_perind_numers.counts or *_perind.counts)')
parser.add_argument('droplist', metavar='droplist', type=str, nargs=1, help='Path to a file containing a list of IDs to exclude from the counts file (typical input is the SamplesToExcludeForPCA.txt file written by /sc/arion/projects/EPIASD/splicingQTL/scripts/pre_fastQTL/PCA/makeRelatedExcludeList.py)')
parser.add_argument('outfile', metavar='outfile', type=str, nargs=1, help='Path to write output file')

args = parser.parse_args()

incountsfile = args.incountsfile[0]
droplist = args.droplist[0]
outfile = args.outfile[0]

incounts = pd.read_csv(incountsfile, sep=' ')
dropme = pd.read_csv(droplist, sep=' ',header=None)

keepme = pd.Series(list(incounts))[~pd.Series(list(incounts)).isin(dropme[0])]

incounts.to_csv(outfile, sep=' ', index=False, columns=keepme)
