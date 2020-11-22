args = commandArgs(trailingOnly = TRUE)

infile = args[1]

print(infile)

qtlres = read.table(infile, head=FALSE, stringsAsFactors=FALSE)

fisher.test(matrix(c(qtlres$V1, qtlres$V2, round(qtlres$V3), qtlres$V2), ncol=2))
