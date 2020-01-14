#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J reformatCovars
#BSUB -P acc_pintod02b
#BSUB -q express
#BSUB -n 8
#BSUB -W 0:10
#BSUB -R rusage[mem=16000]
#BSUB -o /sc/orga/projects/EPIASD/splicingQTL/lsf_output/reformatCovars_plinkPCs_%J_output.txt

ml python/3.7.3

covarsPath='/sc/orga/projects/pintod02b/capstone/leafcutter/data-freeze4/'
covarsFile='leafcutter-input_covar'
metadatFile='/sc/orga/work/belmoj01/sqtl/sampleidmaps/meta_matchedIDs.csv'
plinkFile='/sc/orga/projects/EPIASD/splicingQTL/PCA/plinkPCAresults/Capstone4.sel.idsync.2allele.maf5.mind1.geno1.hwe1e-5.deduped.highLDexcl.indep1500_150_.2/Capstone4.sel.idsync.2allele.maf5.mind1.geno1.hwe1e-5.deduped.highLDexcl.indep1500_150_.2.eigenvec'
relationFile='/sc/orga/projects/EPIASD/splicingQTL/PCA/Capstone4.sel.idsync.2allele.maf5.mind1.geno1.hwe1e-5.relatedness.genome'
dropFile='/sc/orga/projects/EPIASD/splicingQTL/PCA/SamplesToExcludeForPCA.txt'
#dropFile='/sc/orga/projects/EPIASD/splicingQTL/idsync/empty.file'
outPath='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/covar/'

python3 /sc/orga/projects/EPIASD/splicingQTL/idsync/reformatCovars_plinkPCs.py ${covarsPath}${covarsFile}.txt ${metadatFile} ${plinkFile} ${relationFile} ${dropFile} ${outPath}${covarsFile}.20genoPCs.idsync.txt
