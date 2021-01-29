module load ucsc-utils/2015-04-07

chainfile=/hpc/packages/minerva-centos7/liftover/09-Jul-2019/hg19ToHg38.over.chain

indir=$1 # directory containing the bed file to lift over
inbed=$2 # bed file to lift over
outdir=$3 # directory to place the processed files in
addchr=$4 # should 'chr' be prepended to column 1 {true,false}?

mkdir -p ${outdir}

if [ $addchr == 'true' ]
then
  sed 's/^/chr/g' ${indir}${inbed} > ${outdir}${inbed}.hg19withChr
else
  cp ${indir}${inbed} ${outdir}${inbed}.hg19withChr
fi

liftOver -bedPlus=6 ${outdir}${inbed}.hg19withChr ${chainfile} ${outdir}${inbed%.bed}.hg38.bed ${outdir}${inbed%.bed}.unmapped.bed

rm ${outdir}${inbed}.hg19withChr
