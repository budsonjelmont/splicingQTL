#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J reformatCounts
#BSUB -P acc_pintod02b
#BSUB -q premium
#BSUB -n 4
#BSUB -W 04:00
#BSUB -R rusage[mem=16000]
#BSUB -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/reformatCounts_%J_output.txt

ml python/3.7.3

countsdir=/sc/arion/projects/EPIASD/splicingQTL/output/pheno_wasp
countsfile=out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts
metadatfile=/sc/arion/projects/EPIASD/splicingQTL/data/metadata/subsets/meta_matchedIDs_highestRINSampleperSubject_total_1310.tsv
outdir=$countsdir

python3 /sc/arion/projects/EPIASD/splicingQTL/scripts/pre_fastQTL/idsync/reformat_counts.py $countsdir/$countsfile $metadatfile $outdir/${countsfile}.idsync
