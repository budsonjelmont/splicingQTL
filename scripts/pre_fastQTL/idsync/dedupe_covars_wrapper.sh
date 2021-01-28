#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J reformatCovars
#BSUB -P acc_pintod02b
#BSUB -q express
#BSUB -n 8
#BSUB -W 0:30
#BSUB -R rusage[mem=16000]
#BSUB -o /sc/arion/scratch/splicingQTL/lsf_output/reformatCovars_%J_output.txt

ml python/3.7.3

covarsdir=/sc/arion/projects/EPIASD/splicingQTL/output/covars/
covarsfile=leafcutter-input_covar_star_WASP.20genoPCs.idsync.txt
dropfile=/sc/arion/projects/EPIASD/splicingQTL/output/vcf/SamplesToExcludeForPCA.txt
outdir=$covarsdir

python3 dedupe_covars.py $covarsdir/$covarsfile $dropfile $outdir/${covarsfile%.*}.deduped.txt
