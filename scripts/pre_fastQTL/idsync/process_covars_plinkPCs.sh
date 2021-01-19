#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J reformatCovars
#BSUB -P acc_pintod02b
#BSUB -q express
#BSUB -n 8
#BSUB -W 0:10
#BSUB -R rusage[mem=16000]
#BSUB -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/reformatCovars_plinkPCs_%J_output.txt

ml python/3.7.3

covarsPath=/sc/arion/projects/pintod02c/WASP_leafcutter/
covarsFile=leafcutter-input_covar_WASP
metadatFile=/sc/arion/projects/EPIASD/splicingQTL/data/metadata/subsets/meta_matchedIDs_highestRINSampleperSubject_total_1310.tsv
plinkFile=/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/geno_pca_res/Capstone4.sel.hasPhenosOnly.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped--thin-count200000/Capstone4.sel.hasPhenosOnly.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.eigenvec
relationFile=/sc/arion/scratch/belmoj01/splicingQTL/Capstone4.sel.hasPhenosOnly.idsync.2allele.maf01.mind05.geno05.hwe1e-6.relatedness.genome
dropFile=/sc/arion/scratch/belmoj01/splicingQTL/SamplesToExcludeForPCA.txt
outPath=/sc/arion/projects/EPIASD/splicingQTL/output/covar_wasp/

mkdir -p $outPath

python3 /sc/arion/projects/EPIASD/splicingQTL/scripts/pre_fastQTL/idsync/reformat_covars_plinkPCs.py $covarsPath${covarsFile}.txt $metadatFile $plinkFile $relationFile $dropFile $outPath${covarsFile}.20genoPCs.idsync.txt
