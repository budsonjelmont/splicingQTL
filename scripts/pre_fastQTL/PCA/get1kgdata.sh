#!/bin/bash
# To run: bsub [options] < this_scripts_name
#BSUB -J get1kgData
#BSUB -P acc_pintod02b
#BSUB -q normal
#BSUB -n 1
#BSUB -W 8:00
#BSUB -R rusage[mem=70000]
#BSUB -o /sc/orga/projects/EPIASD/splicingQTL/lsf_output/get1KG_data_%J_output.txt

# reference: https://www.biostars.org/p/335605/

cd /sc/orga/projects/EPIASD/splicingQTL/PCA/1kg

# Download the files as VCF.gz (and tab-indices)
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;

for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done

# Download 1000 Genomes PED file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;

# Download the GRCh37 / hg19 reference genome

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
gunzip human_g1k_v37.fasta.gz ;

