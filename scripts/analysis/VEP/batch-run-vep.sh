declare -A invcfs
# Our study
invcfs[OurStudy_sQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/4genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_50HCPs/vcf/all_phase3.dedupeByPos_bestMAF.OurStudy_sQTLSNPs_ALL.vcf
invcfs[OurStudy_allSNPs]=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/geno_wasp/vcf/all_phase3.dedupeByPos_bestMAF.our_allQTLsnps_ALL.vcf
# Walker
invcfs[Walker_sQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/sQTLs/vcf/all_phase3.dedupeByPos_bestMAF.Walker_sQTLSNPs_ALL.vcf
invcfs[Walker_eQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/eQTLs/vcf/all_phase3.dedupeByPos_bestMAF.Walker_eQTLSNPs_ALL.vcf
# Raj
invcfs[Raj_sQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Raj_NatComm_alzh2018/vcf/all_phase3.dedupeByPos_bestMAF.Raj_sQTLSNPs_ALL.vcf
# PEC
invcfs[PEC_eQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECeQTLs/vcf/all_phase3.dedupeByPos_bestMAF.PEC_eQTLSNPs_ALL.vcf
#invcfs[PEC_noneQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECeQTLs/vcf/all_phase3.dedupeByPos_bestMAF.PEC_noneQTLSNPs_ALL.vcf
invcfs[PEC_isoQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoQTLs/hg19/vcf/all_phase3.dedupeByPos_bestMAF.PEC_isoQTLSNPs_ALL.vcf
#invcfs[PEC_nonisoQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECisoQTLs/hg19/vcf/all_phase3.dedupeByPos_bestMAF.PEC_nonisoQTLSNPs_ALL.vcf
invcfs[PEC_tQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECtQTLs/vcf/all_phase3.dedupeByPos_bestMAF.PEC_tQTLSNPs_ALL.vcf
#invcfs[PEC_nontQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PECtQTLs/vcf/all_phase3.dedupeByPos_bestMAF.PEC_nontQTLSNPs_ALL.vcf
invcfs[PEC_allQTLSNPs]=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/PEC/vcf/all_phase3.dedupeByPos_bestMAF.PEC_allQTLsnps_ALL.vcf

resdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/VEP/results_noRegulatory_111620/

mkdir -p $resdir

for dataset in ${!invcfs[@]}
do
 echo "$dataset ..."
 echo "*********************"
 vcf=${invcfs[$dataset]}
 resfile=$dataset
 command="ml vep/97;ml tabix;
 vep -i $vcf \
  -a GRCh37 \
  --hgvs \
  --offline \
  --gene_phenotype \
  --distance 80,80 \
  --most_severe \
  --per gene \
  --cache \
  --dir_cache /sc/arion/work/belmoj01/.vep \
  --fasta /sc/arion/work/belmoj01/.vep/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --force_overwrite \
  --o $resdir$resfile"
 cmd="submitjob 24 -P acc_pintod02b -q premium -c 1 -n 1 -m 60 -J VEP_$dataset ${command}"
 echo $cmd
 $cmd
done
