library(data.table)
library(ggplot2)
library(RColorBrewer)

metadir = '/sc/arion/projects/EPIASD/splicingQTL/data/metadata/'
metafile = 'meta_matchedIDs.csv'
excludefile = '/sc/arion/projects/EPIASD/splicingQTL/scripts/pre_fastQTL/PCA/SamplesToExcludeForPCA.txt' 

meta = fread(paste0(metadir,metafile))
excludes = fread(excludefile, header=F)

# Function to make ancestry plot
vert_group_bar = function(dat,fill,grp,xcol,xname,yname,palette,outpath){
  p=ggplot(data = dat, aes_string(x = xcol, fill = fill), width = 0.6) +
    geom_bar(aes_string(group = grp, fill = fill),colour='black', stat = 'count', position=position_dodge(), width=0.7) +
#    facet_wrap(as.formula(paste0('~',grp))) + # Added as a workaround because group = grp was not working as expected
    scale_fill_brewer('div', palette=palette) +
    xlab(xname) +
    scale_y_continuous(name = yname,expand = c(0,.1)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color='black', size=0.6),
      axis.line.y = element_line(color='black', size=0.6),
      axis.text.x = element_text(colour='black', size=11),
      axis.text.y = element_text(colour='black', size=11),
      axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=12),
      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=12),
      axis.ticks = element_blank(),
      legend.position = c(.75,.875),
      legend.text = element_text(size = 12),
      legend.title=element_blank(),
      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
    )
#  plot(p)
  aspect_ratio = 1.5
  height = 5
  ggsave(outpath,p, height = 7 , width = 7 * aspect_ratio)
}

vert_group_multibar = function(dat,fill,grp,xcol,xname,yname,palette,outpath){
  p=ggplot(data = dat, mapping=aes_string(x = xcol, fill = fill), width = 0.6) +
    #geom_bar(aes_string(group = grp, fill = fill), stat = 'count', position=position_dodge(), width=0.7) +
    geom_bar(stat = 'count', width=0.7) +
    scale_fill_brewer('div', palette=palette) +
    xlab(xname) +
    scale_y_continuous(name = yname,expand = c(0,.1)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(color='black', size=0.6),
      axis.line.y = element_line(color='black', size=0.6),
      axis.text.x = element_text(colour='black', size=15),
      axis.text.y = element_text(colour='black', size=15),
      axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=17),
      axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=17),
      axis.ticks = element_blank(),
      legend.position = c(.75,.875),
      legend.text = element_text(size = 15),
      legend.title=element_blank(),
      plot.margin=unit(c(0.25,1,0.25,1), 'cm')
    )
  aspect_ratio = 1.5
  height = 5
  ggsave(outpath,p, height = 7 , width = 7 * aspect_ratio)
}
vert_multibar = function(df,x,f,name,yaxmax,outdir){
  p = ggplot(data=df, mapping = aes_string(x=x, fill=f)) +
    geom_bar(stat = "count", width = 0.7) +
    scale_fill_manual(values = palette) +
    labs("y" = "Count") +
    ylim(0,yaxmax) + 
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(hjust = 1, angle=80),
      plot.title = element_blank()
    )
  ggsave(paste(outdir,name,'.pdf',sep=''), plot = p, width=5, height=5)
  return(p)
}

# Drop samples that aren't included in fastQTL analysis & drop first column because it is duplicated
meta = meta[!(meta$PSIcountID %in% excludes$V1),]
meta = meta[,-1]

# Set levels (Note: I do this to match the colors in the ancestry PCA plot when the same palette is used in both) 
meta$ethnicity = factor(meta$ethnicity, levels = c('AA','HISP','AS','CAUC','Other'))

# See counts by ancestry & study
table(meta$ethnicity)
table(meta$study)

# Make plots 
# Reported ancestry
vert_group_bar(meta,'ethnicity','ethnicity','ethnicity','Ethnicity','Subjects','Dark2','sQTL_subject_ancestry.png')
vert_group_bar(meta,'ethnicity','ethnicity','ethnicity','Ethnicity','Subjects','Dark2','sQTL_subject_ancestry.pdf')
vert_group_bar(meta,'study','study','study','Study','Subjects','Dark2','sQTL_subjects_per_study.png')
vert_group_bar(meta,'study','study','study','Study','Subjects','Dark2','sQTL_subjects_per_study.pdf')
# Multibar
vert_group_multibar(meta,'ethnicity','study','study','Study','Subjects','Dark2','sQTL_subjects_per_study_mbar_v2.png')
vert_group_multibar(meta,'ethnicity','study','study','Study','Subjects','Dark2','sQTL_subjects_per_study_mbar_v2.pdf')
# Estimated ancestry using 1kg PCA projection
vert_group_bar(meta,'ethnicity_est','ethnicity_est','ethnicity_est','Ethnicity','Subjects','Dark2','sQTL_subject_estimated_ancestry.png')
vert_group_bar(meta,'ethnicity_est','study','study','Study','Subjects','Dark2','sQTL_subject_estimated_ancestry_by_study.png')
