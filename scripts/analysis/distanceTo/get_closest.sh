# Read in a bedfile of sQTL SNPs (fileA) and a bedfile of annotations (fileB) & identify the closest feature in fileB for each SNP in fileA

# To run on sQTL data:
# ./get_closest.sh /sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl.uniqueSNPs.sorted.bed gencode.v19.TSS.PFC_merge_cage_Pitt-Fantom5-119FrontalLob.bed /sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl+ensg.csv nearest_TSS

fileA=$1
fileB=$2
sqtlresfile=$3
newcol=$4

btout="${fileA}_bedtoolsClosest_$(basename ${fileB})"

sqtlresfile_new=${sqtlresfile%.*}+${newcol}.txt

ml bedtools
bedtools closest -D a -a ${fileA} -b ${fileB} | awk -v featname=${newcol} -v distfeatname=${newcol}_dist 'BEGIN{OFS=FS="\t"; print "sid",featname,distfeatname};{print $4,$10,$16}' > ${btout} # -d reports abs val of distance, -D is signed (- value if upstream) 

# Add new columns to sQTL results file 
awk 'BEGIN{OFS=FS="\t"}; NR==FNR{nearestid[$1]=$2; dist[$1]=$3; next};{print $0,nearestid[$6],dist[$6]}' ${btout} ${sqtlresfile} > ${sqtlresfile_new}
