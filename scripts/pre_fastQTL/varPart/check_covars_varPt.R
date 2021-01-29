library(variancePartition)
library(doParallel)
library(preprocessCore)

# Parallelization
cl = makeCluster(25)
registerDoParallel(cl)

# Read input files
countsfile = '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/pheno_wasp/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.deduped.gz.qqnorm_allCombined'
covarsfile = '/sc/arion/projects/EPIASD/splicingQTL/intermediate_files/covar_wasp/leafcutter-input_covar_WASP.20genoPCs.idsync.deduped.10splicingPCs.txt'
dropfile = '/sc/arion/projects/EPIASD/splicingQTL/scripts/pre_fastQTL/varPart/dropTerms.txt'

vpartpath = '/sc/arion/projects/EPIASD/splicingQTL/varPart/'

setwd(vpartpath)

counts = read.table(countsfile, header=TRUE, comment.char='@', stringsAsFactors=FALSE)

# Drop duplicate rows created by cat'ing all the prepare pheno table output
counts=counts[-which(counts$ID=='ID'),]

# Drop first column after setting row names
rownames(counts) = counts$ID
counts = counts[ , !(colnames(counts) %in% colnames(counts)[1:4])]

# Convert all columns to numeric
counts = as.data.frame(sapply(counts, as.numeric))

# Standardize across samples
#counts.sc = apply(counts, 2, scale)

# Normalize across clusters
#counts.sc.norm = normalize.quantiles(counts.sc)

#counts.sc.norm = as.data.frame(counts.sc.norm, row.names=rownames(counts)) 
#colnames(counts.sc.norm) = colnames(counts)

# Write out scaled and normalized data 
#write.table(counts.sc.norm, countsoutfile, row.names=TRUE, col.names=TRUE) 

#counts[] = lapply(counts, convertToNumbers)

covars = read.table(covarsfile, sep='\t', header=TRUE, comment.char='')
rownames(covars) = covars[,1]
covars = covars[,-1]

# Transpose covar data frame
covars = as.data.frame(t(covars))

# Convert non-categorical columns from factors -> numbers
catcols = c('Group','study','sex','individualIDSource','tissue')

for (i in colnames(covars)[!colnames(covars) %in% catcols]){
  covars[,i] = as.numeric(as.character(covars[,i]))
}

covars[,!colnames(covars) %in% catcols] = scale(covars[,!colnames(covars) %in% catcols], center = TRUE, scale = TRUE)  
# Define formula
#form = ~ (1|Group) + ageDeath + ageDeath.squared + (1|study) + (1|sex) + PMI + RIN + RIN.squared + (1|individualIDSource) + (1|tissue) + seqPC3.squared + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC15 + seqPC16 + seqPC17 + seqPC18 + seqPC19 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC26 + seqPC27 + seqPC28 + seqPC29 + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + HCP1 + HCP2 + HCP3 + HCP4 + HCP5 + HCP6 + HCP7 + HCP8 + HCP9 + HCP10 + HCP11 + HCP12 + HCP13 + HCP14 + HCP15 + HCP16 + HCP17 + HCP18 + HCP19 + HCP20 # ALL covariates
form = ~ (1|Group) + ageDeath + ageDeath.squared + (1|study) + (1|sex) + PMI + RIN + RIN.squared + (1|individualIDSource) + (1|tissue) + seqPC3.squared + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC15 + seqPC16 + seqPC17 + seqPC18 + seqPC19 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC26 + seqPC27 + seqPC28 + seqPC29 + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + genotypePC6 + genotypePC7 + genotypePC8 + genotypePC9 + genotypePC10 + genotypePC11 + genotypePC12 + genotypePC13 + genotypePC14 + genotypePC15 + genotypePC16 + genotypePC17 + genotypePC18 + genotypePC19 + genotypePC20 + splicingPC1 + splicingPC2 + splicingPC3 + splicingPC4 + splicingPC5 + splicingPC6 + splicingPC7 + splicingPC8 + splicingPC9 + splicingPC10 # ALL covariates, no HCPs
#form = ~ (1|Group) + ageDeath + (1|sex) + PMI + RIN.squared + (1|individualIDSource) + (1|tissue) + seqPC3.squared + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC15 + seqPC16 + seqPC17 + seqPC18 + seqPC19 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC26 + seqPC27 + seqPC28 + seqPC29 + genotypePC1 + genotypePC2 + HCP1 + HCP2 + HCP3 + HCP4 + HCP5 + HCP6 + HCP7 + HCP8 + HCP9 + HCP10 # ALL covariates except correlated
#form = ~ (1|Group) + ageDeath + (1|study) + (1|sex) + PMI + RIN + (1|tissue) + seqPC3.squared + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC15 + seqPC16 + seqPC17 + seqPC18 + seqPC19 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC26 + seqPC27 + seqPC28 + seqPC29 + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + splicingPC3 + splicingPC6 + splicingPC7 + splicingPC8 + splicingPC9 + splicingPC10 + HCP1 + HCP2 + HCP3 + HCP4 + HCP5 + HCP6 + HCP7 + HCP8 + HCP9 + HCP10 + HCP11 + HCP12 + HCP13 + HCP14 + HCP15 + HCP16 + HCP17 + HCP18 + HCP19 + HCP20

#form = update.formula(form, )

# Fit model and extract results
varPart = fitExtractVarPartModel( counts, form, covars ) 

# Save workspace image  
save.image()

# Sort variables by median fraction of variance explained
vp = sortCols( varPart )
# Make data frame for easy access to variance partition stats
vp_summary = data.frame(do.call(cbind, lapply(vp, summary)))

# Figure 1a# Bar plot of variance fractions for the first 20 genes
fig = plotPercentBars( vp[1:20,] )
ggsave(paste(vpartpath,'varpart_barplot.png',sep=''), fig)
## Figure 1b# violin plot of contribution of each variable to total variance
vp_plot = vp
colnames(vp_plot) = paste(colnames(vp),' (',as.character(signif(vp_summary['Median',],3)*100),'%)',sep='')
fig = plotVarPart( vp_plot )
#fig = plotVarPart( vp_plot , xlab=paste(colnames(vp),' (',as.character(vp['Median',]),'%)',sep=''))
fig = fig + theme(axis.text.x = element_text(size = 7, angle = 40))
ggsave(paste(vpartpath,'varpart_vplot.png',sep=''), fig)

# Reformat formula to exclude random effects (since canCorPairs currently can't handle them)
#reform = ~ Group + ageDeath + ageDeath.squared + study + sex + PMI + RIN + RIN.squared + individualIDSource + tissue + seqPC3.squared + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC15 + seqPC16 + seqPC17 + seqPC18 + seqPC19 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC26 + seqPC27 + seqPC28 + seqPC29 + SV1 + SV2 + SV3 + SV4 + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + genotypePC6 + genotypePC7 + genotypePC8 + genotypePC9 + genotypePC10 + genotypePC11 + genotypePC12 + genotypePC13 + genotypePC14 + genotypePC15 + genotypePC16 + genotypePC17 + genotypePC18 + genotypePC19 + genotypePC20 + splicingPC1 + splicingPC2 + splicingPC3 + splicingPC4 + splicingPC5 + splicingPC6 + splicingPC7 + splicingPC8 + splicingPC9 + splicingPC10 + HCP1 + HCP2 + HCP3 + HCP4 + HCP5 + HCP6 + HCP7 + HCP8 + HCP9 + HCP10 + HCP11 + HCP12 + HCP13 + HCP14 + HCP15 + HCP16 + HCP17 + HCP18 + HCP19 + HCP20
reform = ~ Group + ageDeath + ageDeath.squared + study + sex + PMI + RIN + RIN.squared + individualIDSource + tissue + seqPC3.squared + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC15 + seqPC16 + seqPC17 + seqPC18 + seqPC19 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC26 + seqPC27 + seqPC28 + seqPC29 + SV1 + SV2 + SV3 + SV4 + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + genotypePC6 + genotypePC7 + genotypePC8 + genotypePC9 + genotypePC10 + genotypePC11 + genotypePC12 + genotypePC13 + genotypePC14 + genotypePC15 + genotypePC16 + genotypePC17 + genotypePC18 + genotypePC19 + genotypePC20 + splicingPC1 + splicingPC2 + splicingPC3 + splicingPC4 + splicingPC5 + splicingPC6 + splicingPC7 + splicingPC8 + splicingPC9 + splicingPC10 # All covars, no HCPs 

# Assess correlation between covariates
covarCor = canCorPairs(reform, covars, showWarnings = TRUE)

png(paste(vpartpath,'varpart_covarcorplot_hclust.png',sep=''),
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4)
plotCorrMatrix(covarCor, margins= c(16.5,8), dendrogram='none',sort=F)
dev.off()

# Plot correlation between covars
library(reshape2)
melted_covarCor =  melt(covarCor)

#reorder_cormat = function(cormat){
#  # Use correlation between variables as distance
#  dd = as.dist((1-cormat)/2)
#  hc = hclust(dd)
#  cormat = cormat[hc$order, hc$order]
#}

corfig = ggplot(data = melted_covarCor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "white", high = "red", mid="pink",
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
   size = 8, hjust = 1))+
  coord_fixed()

ggsave(paste(vpartpath,'varpart_covarcorplot.png',sep=''), corfig)

# Make list of covariates to drop (correlation > X and explains less variability)
x = 0.9
covarcols = colnames(vp)[1:length(colnames(vp))-1]
todrop = lapply(covarcols, function(covar1){
  lapply(covarcols, function(covar2,covar1){
    print(paste(covar1,covar2,sep='    '))
    if (covar1 == covar2){return(NA)}
    else {
      tryCatch(
        {if (covarCor[covar1,covar2] > x){
          if(median(vp[,covar1]) > median(vp[,covar2])){
            return(covar2)
          } else if (median(vp[,covar1]) < median(vp[,covar2])){
            return(covar1)
          } else {return('???')}
        } else {return(NA)}
        },
        error = function(c){return(NA)}
      )
    }
  },covar1)
  }
)

todrop = unlist(unlist(todrop))
todrop = unique(todrop[!is.na(todrop)])

write.table(todrop,paste('covarstodrop_corrthresh',as.character(x),'.txt',sep=''),row.names=F,quote=F,col.names=F)

# Write the data table out
write.table(varPart, file = paste(vpartpath,'varPart',sep=''), append = FALSE, quote = FALSE, sep = '\t', eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
