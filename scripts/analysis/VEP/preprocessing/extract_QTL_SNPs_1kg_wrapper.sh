#!/bin/bash

declare -A DATASET # 1st dimension array. One per dataset.
declare -A DATASETS # 2nd dimension array. Each key is of format [integer_id_of_1st_d_array],[array_key_from_ast_d_array]
i=1 # Dataset counter

# Call after each DATASET declaration to populate 2d-array
function populate_datasets(){
  for key in "${!DATASET[@]}"; do
    DATASETS[$i,$key]=${DATASET[$key]}
  done
  ((i++))
}

# All our SNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/
DATASET["INSNPSFILE"]=Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY.snplist
DATASET["DATATAG"]=OurStudy_allSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# Our sSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/output/sqtls/minCovars+seqPC9/4genoPCs/35HCPs/qtls/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=OurStudy_sQTLSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# Walker sSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/Walker2018/sQTLs/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=our_sQTLSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# Walker eSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/Walker2018/eQTLs/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=Walker_eQTLSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# Raj sSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/Raj2018/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=Raj_sQTLSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# PEC all SNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/PEC/
DATASET["INSNPSFILE"]=SNP_Information_Table_with_Alleles.chrpos.txt
DATASET["DATATAG"]=PEC_allSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# PEC isoSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/PEC/PECisoQTLs/hg19/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=PEC_isoSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# PEC eSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/PEC/PECeQTLs/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=PEC_eSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# PEC tSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/PEC/PECtQTLs/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=PEC_tSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# PEC cSNPs
DATASET["INSNPSDIR"]=/sc/arion/projects/EPIASD/splicingQTL/ext_data/QTLs/PEC/PECcQTLs/
DATASET["INSNPSFILE"]=qtl_snps_unique.txt
DATASET["DATATAG"]=PEC_cSNPs_ALL
DATASET["REFDIR"]=/sc/arion/scratch/belmoj01/QTL_VCF/
DATASET["TMPDIR"]=${DATASET["INSNPSDIR"]}vcf/
DATASET["REFNAME"]=all_phase3.dedupeByPos_bestMAF
populate_datasets 

# Loop and call SNP extraction script 
for j in $(seq 1 $((i-1)))
do
 echo "################################## Processing ${DATASETS[$j,"DATATAG"]} ##################################"
 echo "./extract_QTL_SNPs_1kg.sh ${DATASETS[$j,"INSNPSDIR"]} ${DATASETS[$j,"INSNPSFILE"]} ${DATASETS[$j,"DATATAG"]} ${DATASETS[$j,"REFDIR"]} ${DATASETS[$j,"TMPDIR"]} ${DATASETS[$j,"REFNAME"]}"
 ./extract_QTL_SNPs_1kg.sh ${DATASETS[$j,"INSNPSDIR"]} ${DATASETS[$j,"INSNPSFILE"]} ${DATASETS[$j,"DATATAG"]} ${DATASETS[$j,"REFDIR"]} ${DATASETS[$j,"TMPDIR"]} ${DATASETS[$j,"REFNAME"]}
done
