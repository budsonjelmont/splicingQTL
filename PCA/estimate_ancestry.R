library(plinkQC)
#library(RColorBrewer)

setwd('~/Documents/sqtl/PCA/')

extdatdir = system.file("extdata", package="plinkQC")
indir = '~/Documents/sqtl/PCA/ASD-EPI_Plates1-2'
name = 'ASD-EPI_Plates1-2'
refname = '1kg_phase3'
refSamplesFile= '~/Documents/sqtl/PCA/20130606_g1k.ped_refindividpop.txt'
prefixMergedDataset = paste(name, ".", refname, sep="")

# Read in reference samples file & assign colors
refsamples = read.table(refSamplesFile, sep='\t', header=T)
refcolors = data.frame(Pop=unique(refsamples['Pop']))
refcolors['Color'] = palette(rainbow(nrow(refcolors)))

# Return value: Named [list] with 
#  i) fail_ancestry, containing a [data.frame] with FID and IID of non-European individuals and 
#  ii) p_ancestry, a ggplot2-object 'containing' a scatter plot of PC1 versus PC2 colour-coded for samples of the reference populations and the study population.
exclude_ancestry =
  evaluate_check_ancestry(indir=indir, 
    name=name,
    prefixMergedDataset=prefixMergedDataset,
    refSamplesFile=refSamplesFile,
  #                          refColorsFile=paste(extdatdir, "/HapMap_PopColors.txt",
  #                                             sep=""),
    refColors=refcolors,
    studyColor='black',
    interactive=TRUE
  )

# Save plot
library(ggplot2)
pwidth = 13 # Width of plot (may need to increase if # of distinct pops in graph is high)
ggsave(filename = paste0(indir,'/',name,'_merge_1kg_PCA.png'),
         plot = exclude_ancestry$p_ancestry,
         device = 'png',
         width = pwidth
       )
# Spit out list of non-European individuals
write.table(exclude_ancestry$fail_ancestry,
             paste0(indir,'/',name,'_merge_1kg_PCA-nonEuroIndivs.txt'),
             quote=F,
             sep='\t'
            )
