# reference: https://www.biostars.org/p/335605/

cd /sc/orga/projects/EPIASD/splicingQTL/PCA/1kg

# Convert the BCF files to PLINK format
for chr in {1..22}; do
    plink --noweb --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed \
    --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done

# Prune variants from each chromosome

mkdir Pruned ;

for chr in {1..22}; do
    plink --noweb --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --maf 0.10 --indep 1500 50 0.2 \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;

    plink --noweb --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.prune.in --make-bed \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done
