# Read in the 1kg VCF already filtered to include only study SNPs & remap the pos to the lifted over hg38 pos. Leave the chr:pos ID the same so that we can match back to the original sQTL results.
ml plink/1.90
ml bcftools

vcfindir=/sc/arion/scratch/belmoj01/QTL_VCF/
vcfin=all_phase3.dedupeByPos_bestMAF.OurStudy_sQTLSNPs_ALL.vcf
qtlshg38in=/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output_wasp/4genoPCs_nogenoInHCP/deduped_mincovars+seqPC9_50HCPs/qtls+hg38.txt
vcfoutdir=$vcfindir
vcfout=${vcfin%.*}.hg38pos

idhg38posmap=${qtlshg38in%.*}.id_hg38pos.map
idsliftedover=${qtlshg38in%.*}.id_liftedover.tsv

awk -v sidcol=sid -v hg38poscol=sid_end_hg38 'BEGIN{OFS=FS="\t"};
   (NR==1){sidcolnum=-1; hg38poscolnum=-1;
    for(i=1;i<=NF;i++){
      if($(i)==sidcol){sidcolnum=i};
      if($(i)==hg38poscol){hg38poscolnum=i};
    };
   } 
   (NR!=1){
    if($hg38poscolnum!=""){
      print $sidcolnum,$hg38poscolnum 
    };
   }
  ' $qtlshg38in > $idhg38posmap.withdups

# Remove duplicates SNP IDs from the file produced in the previous step
sort $idhg38posmap.withdups | uniq > $idhg38posmap

# Read in VCF & throw out any SNPs that didn't lift over to hg38. Then remap positions to hg38 coordinates (but don't change the chr:pos IDs).

plink --vcf $vcfindir$vcfin \
  --double-id \
  --update-map $idhg38posmap \
  --extract $idhg38posmap \
  --make-bed \
  --out $vcfoutdir$vcfout.unsorted

# Call plink one more time to sort IDs before zipping & indexing
plink --bfile $vcfoutdir$vcfout.unsorted \
  --recode vcf-iid \
  --out $vcfoutdir$vcfout

# Zip & tabix
bgzip $vcfoutdir$vcfout.vcf
tabix -p vcf $vcfoutdir$vcfout.vcf.gz

# Remove temp files we made
rm $idhg38posmap.withdups
rm $idhg38posmap
