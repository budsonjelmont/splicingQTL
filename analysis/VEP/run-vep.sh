#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J VEP
#BSUB -P acc_pintod02b
#BSUB -q express
#BSUB -n 1
#BSUB -W 2:00
#BSUB -R rusage[mem=100000]
#BSUB -o /sc/arion/scratch/belmoj01/VEPannotation_%J_output.txt

vcfdir=/sc/arion/scratch/belmoj01/QTL_VCF/
vcffile=all_phase3.dedupeByPos_bestMAF.our_nonQTLSNPs.vcf
resfile=our_nonQTLSNPs_offlineRun
resdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results/

ml vep/97
ml tabix

mkdir -p $resdir

vep -i $vcfdir$vcffile \
  -a GRCh37 \
  --hgvs \
  --offline \
  --gene_phenotype \
  --regulatory \
  --distance 80,80 \
  --most_severe \
  --per gene \
  --cache \
  --dir_cache /sc/hydra/work/belmoj01/.vep/ \
  --fasta /sc/hydra/work/belmoj01/.vep/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --force_overwrite \
  --o $resdir$resfile
