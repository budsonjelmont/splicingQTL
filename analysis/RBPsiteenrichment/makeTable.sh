# Read in all .txt files in the RBP enrichment results directory and compile them into a single table 
resdir='/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/hg38_liftOver/qtlTools_RBP_results/'
restable=RBPenrichmentTable.tsv

resfiles=(${resdir}*.txt);

rm -f ${resdir}${restable}

# Split filenames to get RBP tag & add it to each line of results table
for res in ${resfiles[@]}
do
  rbp=$(echo ${res} | cut -d'-' -f 2 | cut -d'.' -f 1)
  echo ${res};
  echo ${rbp};
  awk -v rbp=${rbp} 'BEGIN{OFS=FS="\t"};{print rbp,$0};' ${res} >> ${resdir}${restable}
done
