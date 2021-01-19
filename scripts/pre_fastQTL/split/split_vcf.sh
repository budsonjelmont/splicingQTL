idlist=$1
idcol=$2
invcf=$3
outdir=$4

mkdir -p $outdir
outbase=$(basename $invcf)

echo $outbase

while read line
do
#  echo $line
  arr=($(echo $line | tr ' ' "\n"))
  echo "${arr[$idcol]} ${arr[$idcol]}" > $outdir/${arr[$idcol]}.txt
  cmd="ml plink/1.90; plink --vcf $invcf --double-id --keep $outdir/${arr[$idcol]}.txt --recode vcf-iid --out $outdir/${outbase}.${arr[$idcol]}"
  echo $cmd
  submitjob 1 -P acc_pintod02b -q express -c 1 -n 1 -m 30 -J splitVCF_${arr[$idcol]} $cmd 
done < $idlist
