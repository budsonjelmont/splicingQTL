#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J plinkRemoveLinkedSNPs
#BSUB -P acc_pintod02b
#BSUB -q express
#BSUB -n 1
#BSUB -W 00:45
#BSUB -R rusage[mem=50000]
#BSUB -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/remove_LD_SNPs_%J_output.txt

# Preprocess 1KG reference to make a EUR-only reference containing only SNPs with no duplicate RSIDs/chr:pos that
# can then be used to calculate LD & extract SNPs for GREGOR 

module load plink/1.90

inSNPsdir=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/4genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_50HCPs/
#inSNPsdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/sQTLs/
#inSNPsdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Raj_NatComm_alzh2018/
#inSNPsdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/eQTLs/
#inSNPsdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoQTLs/hg19/

# File listing chr(number only):pos for all SNPs to extract from reference
#inSNPsfile=uniqueLeadSNPuniquephenos-chrpos.txt
inSNPsfile=uniquenonsigSNPs-chrpos.txt
#inSNPsfile=uniqueSNPs-chrpos.txt

dataset=OurStudy_nonQTLSNPs_ALL
refdir=/sc/arion/scratch/belmoj01/QTL_VCF/
tmpdir=$refdir
refname=all_phase3.dedupeByPos_bestMAF
outname=$refname.$dataset

##################################################

# Make directory to do our dirty work in
mkdir -p $tmpdir

##################################################

# Make 3 col BED file of coordinates 
sort $inSNPsdir$inSNPsfile | uniq | gawk 'BEGIN{FS="\t"}; {match($1,"([0-9XYM]+):([0-9]+)", coords); print coords[1],coords[2],coords[2],coords[1]":"coords[2]}' > $inSNPsdir$inSNPsfile.extract.in

# Extract dataset SNPs & recode the resulting VCF   
plink --bfile $refdir$refname \
      --allow-extra-chr \
      --extract range $inSNPsdir$inSNPsfile.extract.in \
      --recode vcf \
      --out $tmpdir$outname 

# Clean up after ourselves
#mv plink.log ${inSNPsdir}$filetag_plink.log
#rm plink*
