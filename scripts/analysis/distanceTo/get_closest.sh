# Read in a bedfile of sQTL SNPs (fileA) and a bedfile of annotations (fileB) & identify the closest feature in fileB for each SNP in fileA

# To run on sQTL data:
# ./get_closest.sh /sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/sqtlSNPS_unique.sorted.bed /sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/TSS/gencode.v19.annotation_capped_sites_nr_with_confidence.reformat.PFConly.sorted.bed /sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/20genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_15HCPs/significant_sqtl+ensg.txt nearest_TSS sid

fileA=$1 # 6-column bed file
fileB=$2 # 6 or 9-column bed file 
sqtlresfile=$3 # sQTL results file to add new column to
newcol=$4 # name of new column
resfileidcol=$5 # name of the ID column in the results file (e.g. sid for fastQTL output)

btout="${fileA}_bedtoolsClosest_$(basename ${fileB})"

sqtlresfile_new=${sqtlresfile%.*}+${newcol}.txt

# Return the number of a column in tab-delimited file given its header label
function get_col_num(){
 file=$1
 colname=$2
 sed -n $'1s/\t/\\\n/gp' ${file} | grep -nx ${colname} | cut -d: -f1
}

module load bedtools
bedtools closest -D a -a ${fileA} -b ${fileB} | awk -v idcol=${resfileidcol} -v featname=${newcol} -v distfeatname=${newcol}_dist 'BEGIN{OFS=FS="\t"; print idcol,featname,distfeatname};{print $4,$10,$NF}' > ${btout} # -d reports abs val of distance, -D is signed (- value if upstream) 

# Add new columns to sQTL results file 
awk -v idcol=$(get_col_num ${sqtlresfile} ${resfileidcol}) 'BEGIN{OFS=FS="\t"}; NR==FNR{nearestid[$1]=$2; dist[$1]=$3; next};{print $0,nearestid[$idcol],dist[$idcol]}' ${btout} ${sqtlresfile} > ${sqtlresfile_new}
