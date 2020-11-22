classfile=/sc/hydra/projects/pintod02c/isoseq_MAP/data_freeze_070920/map-paper-merge-capped-a100-cagedistfix_classification_extra_idfixed_suppa1vsAll-PFCoutHSBin.txt
column=AF_FBvsHSB,CBL,CTX
#column=AF_FBvsHSB\,CBL\,CTX
pattern=^chr

awk -v targetcol=${column} -v targetpattern=${pattern} '(NR==1){targetcolnum=-1;
 for(i=1;i<=NF;i++)
  if($(i)==targetcol){targetcolnum=i};
   next
  };
  {print $targetcolnum};
' $classfile > out.txt
#  $targetcolnum ~ /$targetpattern/ {print $0};' $classfile > $(basename ${classfile} .txt)_${column}Selected.txt
