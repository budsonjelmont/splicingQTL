# Read in a pre-filtered classification file (or any file with #isoform in the header) and extract all PBIDs. Then extract those PBIDs from the full IsoSeq map bedfile

metafile=$1 # Metadata file to extract PBIDs from (look for column '#isoform') that contains all and only the PBIDs of interest # e.g. /sc/hydra/work/belmoj01/lncrna_map/isoseq_map_analysis/map_summaries/PFC-merge-freeze_stats-070920/anySource/codingandnon/monoandmulti/transcripts_lncRNA.tsv
inbed=$2 # Bedfile containing PBIDs in the 4th column. The rows matching PBIDs in the #isoform column of the classification file will be extracted e.g. /sc/hydra/projects/pintod02c/isoseq_MAP/track_hub2/tama_capped_a100_cagedistfix/map-paper-merge-capped-a100-cagedistfix_corrected.bed # BEDfile containing all transcripts from our map
outbed=$3 # Path & name of output bedfile

gawk -v pbidcol='#isoform' 'BEGIN{OFS=FS="\t"};(NR==1 && NR==FNR){pbidcolnum=-1;
    for(i=1;i<=NF;i++)
      if($(i)==pbidcol){pbidcolnum=i};
      next
    };(NR!=1 && NR==FNR){PBarr[$pbidcolnum];next};
    (NR!=FNR){match($4,"(PB[^|]*)", pbid); if(pbid[1] in PBarr)print $0;}' ${metafile} ${inbed} > ${outbed}.tmp

# Remove transcripts mapping to non-canonical chromosomes
grep -E '^chr[0-9XY]{1,2}' ${outbed}.tmp > ${outbed}

# Remove temp file
rm ${outbed}.tmp
