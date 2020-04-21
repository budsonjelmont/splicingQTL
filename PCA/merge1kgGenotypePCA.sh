# Source: https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html
ml plink

datdir='/sc/arion/projects/EPIASD/splicingQTL/PCA/studydata'
name='ASD-EPI_Plates1-2'
refdir='/sc/arion/projects/EPIASD/splicingQTL/PCA/1kg_phase3'
refname='all_phase3'
reference='1kg_phase3'
highld='high-LD-regions-hg38-GRCh38.txt'

mkdir -r $datdir/plink_log

# Prune study data by pruning sites in LD & also removing pre-computed high-LD areas
plink --bfile  $datdir/$name \
      --exclude range  $refdir/$highld \
      --indep-pairwise 50 5 0.2 \
      --out $datdir/$name
mv  $datdir/$name.prune.log $datdir/plink_log/$name.prune

plink --bfile  $datdir/$name \
      --extract $datdir/$name.prune.in \
      --make-bed \
      --out $datdir/$name.pruned
mv  $datdir/$name.pruned.log $datdir/plink_log/$name.pruned

# Filter reference data for the same SNP set as in study
plink --bfile  $refdir/$refname \
      --allow-extra-chr \
      --extract $datdir/$name.prune.in \
      --make-bed \
      --out $datdir/$refname.pruned
mv  $datdir/$refname.pruned.log $datdir/plink_log/$refname.pruned

# Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $datdir/$name.pruned.bim $datdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $datdir/$refname.toUpdateChr

plink --bfile $datdir/$refname.pruned \
      --update-chr $datdir/$refname.toUpdateChr 1 2 \
      --make-bed \
      --out $datdir/$refname.updateChr
mv $datdir/$refname.updateChr.log $datdir/plink_log/$refname.updateChr.log

# Find variants with position mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $datdir/$name.pruned.bim $datdir/$refname.pruned.bim > \
    $datdir/${refname}.toUpdatePos

# Find allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $datdir/$name.pruned.bim $datdir/$refname.pruned.bim > \
    $datdir/$refname.toFlip

# Update positions and flip alleles
plink --bfile $datdir/$refname.updateChr \
      --update-map $datdir/$refname.toUpdatePos 1 2 \
      --flip $datdir/$refname.toFlip \
      --make-bed \
      --out $datdir/$refname.flipped
mv $datdir/$refname.flipped.log $datdir/plink_log/$refname.flipped.log

# Remove mismatches
# Any alleles that do not match after allele flipping, are identified and removed from the reference dataset.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $datdir/$name.pruned.bim $datdir/$refname.flipped.bim > \
    $datdir/$refname.mismatch

plink --bfile $datdir/$refname.flipped \
      --exclude $datdir/$refname.mismatch \
      --make-bed \
      --out $datdir/$refname.clean
mv $datdir/$refname.clean.log $datdir/plink_log/$refname.clean.log

# Merge study genotypes and reference data
# The matching study and reference dataset can now be merged into a combined dataset with plink –bmerge. If all steps outlined above were conducted successfully, no mismatch errors should occur.

plink --bfile $datdir/$name.pruned  \
      --bmerge $datdir/$refname.clean.bed $datdir/$refname.clean.bim \
         $datdir/$refname.clean.fam  \
      --make-bed \
      --out $datdir/$name.merge.$refname
mv $datdir/$name.merge.$refname.log $datdir/plink_log

# PCA on the merged data
# We can now run principal component analysis on the combined dataset using plink –pca which returns a .eigenvec file with the family and individual ID in columns 1 and 2, followed by the first 20 principal components.

plink --bfile $datdir/$name.merge.$refname \
      --pca \
      --out $datdir/$name.$reference
mv $datdir/$name.$reference.log $datdir/plink_log
