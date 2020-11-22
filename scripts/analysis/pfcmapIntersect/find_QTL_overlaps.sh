ml bedtools

qtldir=/sc/hydra/projects/pintod02c/hadjie01_wd/QTLs/all_QTLs
qtlfilesbase=(our_sQTLs Walker_sQTLs Walker_eQTLs PEC_eQTLs PEC_isoQTLs PEC_tQTLs PEC_cQTLs)

gwasdir=/sc/hydra/projects/pintod02c/hadjie01_wd/QTLs/all_GWAS
gwasfilesbase=(Grove2019_GWAS ILAE2018_GWAS PGC2020_GWAS)

mapdir=/sc/hydra/projects/pintod02c/hadjie01_wd/QTLs/all_maps
mapfilesbase=(CBL_known CBL_novel PFC_known PFC_novel)

resdir=/sc/hydra/projects/pintod02c/hadjie01_wd/QTLs/map_Overlaps
statsfile=$resdir/overlap_stats.tsv

mkdir -p $resdir

rm $statsfile

echo 'comparison	n_overlap_transcripts	n_overlap_genes	n_total_transcripts	n_total_genes' > $statsfile

for map in ${mapfilesbase[@]}
do
  echo "$map ..."
  for qtl in ${qtlfilesbase[@]}
  do
  echo "  $qtl ..."
   bedtools intersect -wo -a $mapdir/$map.bed -b $qtldir/$qtl.bed > $resdir/${map}_overlaps_$qtl
   python3 report_unique.py $resdir/${map}_overlaps_$qtl $mapdir/$map.bed $qtldir/$qtl.bed $statsfile ${map}_overlaps_$qtl 
  done
  for gwas in ${gwasfilesbase[@]}
  do
  echo "  $gwas ..."
   bedtools intersect -wo -a $mapdir/$map.bed -b $gwasdir/$gwas.bed > $resdir/${map}_overlaps_$gwas
   python3 report_unique.py $resdir/${map}_overlaps_$gwas $mapdir/$map.bed $gwasdir/$gwas.bed $statsfile ${map}_overlaps_$gwas 
  done
done
