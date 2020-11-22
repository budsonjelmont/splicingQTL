# SE for QTL data is not returned by fastQTL, so we must calculate it. See discussion here: https://github.com/stephenslab/gtexresults/issues/5

library(coloc)

basedir = '/sc/hydra/projects/pintod02c/datasets-external/disease-lists/PGC2020/loci_to_test/'
dirpattern = 'loci_\\.*'
dirs = list.files(path = basedir, pattern = dirpattern, recursive = F, all.files=F, include.dirs=T, full.names=F)

#  Convert data frame to named list by column
listify = function(DF,type){
  if(type=='quant'){
    list(snp=DF$snp, pvalues=DF$pvalues, beta=DF$beta, varbeta=DF$varbeta, MAF=DF$MAF, N=DF$N[1],type='quant')
  } else if (type=='cc'){
    list(snp=DF$snp, pvalues=DF$pvalues, beta=DF$beta, varbeta=DF$varbeta, MAF=DF$MAF, N=DF$N[1], s=DF$s[1],type='cc')
  }
}

# Loop over directories, and if there are files in them run coloc 
coloc.df.list = lapply(dirs,function(dir){
  gwasfile = paste0(basedir,dir,'/GWAS_SNPs.coloc.dat')
  qtlfile = paste0(basedir,dir,'/QTL_SNPs.coloc.dat')
  empty.df = t(data.frame(c('nsnps'=0, 'PP.H0.abf'=NA, 'PP.H1.abf'=NA, 'PP.H2.abf'=NA, 'PP.H3.abf'=NA, 'PP.H4.abf'=NA))) # Returned if datasets have no SNPs in common/no QTL SNPs are found in this locus
  rownames(empty.df) = dir
  if (file.exists(gwasfile) & file.exists(gwasfile)){
    gwas = read.table(gwasfile, sep='\t', header=T, as.is=T) 
    qtl = read.table(qtlfile, sep='\t', header=T, as.is=T) 
    dat1 = qtl
    dat1$type = 'quant'
    dat2 = gwas
    dat2$type = 'cc' 
    # Check for any overlapping SNPs before attempting coloc
    if(any(dat1$snp %in% dat2$snp)){
      dataset1=listify(dat1,'quant')
      dataset2=listify(dat2,'cc')
      dataset1$varbeta[dataset1$varbeta==0]=0.0000001 # This adjustment is necessary because of records where varbeta = 0, which throws errors in coloc. This may be due to a precision error in the QTL results processing.
      coloc.res = coloc.abf(dataset1, dataset2, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
      write.table(coloc.res$results, paste0(basedir,dir,'/coloc.dat'), sep='\t', quote=F, row.names=F, col.names=T)
      summary = t(as.data.frame(coloc.res$summary))
      rownames(summary) = dir
      return(summary)
    } else {  
     return(empty.df)
    }
  } else {
    return(empty.df)
  }
})

coloc.df = do.call(rbind, coloc.df.list)
write.table(coloc.df, paste0(basedir,'coloc.res.dat'),sep='\t',quote=F,row.names=T,col.names=T)
write.table(coloc.df[which(coloc.df[,'PP.H4.abf'] > 0.95),], paste0(basedir,'coloc.res_PPH4gt95.dat'),sep='\t',quote=F,row.names=T,col.names=T)
