fenrichPath='/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/fqtl_output/20genoPCs_nogeno/deduped_mincovars_40HCPs/'
fenrichOut='RBP_qtltools_enrich.txt'

D=read.table(paste(fenrichPath,fenrichOut,sep=''), head=FALSE, stringsAsFactors=FALSE)
fisher.test(matrix(c(D$V1, D$V2, round(D$V3), D$V2), ncol=2))
