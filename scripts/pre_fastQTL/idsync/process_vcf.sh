# Extract subject IDs (1310 total unique) from the Capstone4 vcf who have RNAseq data

#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J reformatvcf
#BSUB -P acc_pintod02b
#BSUB -q express
#BSUB -n 1
#BSUB -W 00:50
#BSUB -R rusage[mem=200000]
#BSUB -o /sc/arion/scratch/belmoj01/reformatvcf_%J_output.txt

module load python/3.7.3
module load plink/1.90 

oldvcfpath=/sc/arion/projects/EPIASD/PsychENCODE_phase12_genotypes/Capstone4VCFreformat/
oldvcffile=Capstone4.sel
metadatafile=/sc/arion/projects/EPIASD/splicingQTL/metadata/meta_matchedIDs.csv
seqedvcffile=${oldvcffile}.hasPhenosOnly
newvcfpath=/sc/arion/scratch/belmoj01/splicingQTL/
newvcffile=${seqedvcffile}.idsync
extractme=/sc/arion/projects/EPIASD/splicingQTL/data/metadata/subsets/meta_matchedIDs_highestRINSampleperSubject_total_1310.idmap

mkdir -p $newvcfpath

# Extract subjects who have RNAseq data 
plink --vcf $oldvcfpath/${oldvcffile}.vcf --double-id --keep $extractme --recode vcf-iid --out $newvcfpath/$seqedvcffile
# Remap their IDs to be consistent with those found in the leafcutter counts file
plink --vcf $newvcfpath/$seqedvcffile.vcf --double-id --update-ids $extractme --recode vcf-iid --out $newvcfpath/$newvcffile

