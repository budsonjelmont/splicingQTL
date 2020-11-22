# Read in a list of RBPs & search for a BED file from POSTAR containing the binding sites for each

rbpstested=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Raj_NatComm_alzh2018/RBPsTested.txt
beddir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/humanRBPsites/
filebase=human_RBP_binding_sites-
gregorbedfilelist=/sc/arion/projects/EPIASD/splicingQTL/analysis/GREGOR/bedindexes/RajCLIPDBRBPs.bed.list.txt

IFS=$'\n' read -d '' -r -a rbps < $rbpstested 

rm $gregorbedfilelist

for r in ${rbps[@]}
do
  echo ${r};
  if test -f "${beddir}${filebase}${r}.bed"; then
    ls ${beddir}${filebase}${r}.bed >> $gregorbedfilelist
  fi
done
