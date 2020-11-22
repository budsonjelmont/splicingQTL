options(scipen=100, digits=3)

# read in the eigenvectors, produced in PLINK
eigenvec <- data.frame(read.table("plink.eigenvec", header=FALSE, skip=0, sep=" "))
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste("Principal Component ", c(1:20), sep="")

# read in the PED data
PED <- data.frame(read.table("20130606_g1k.ped", header=TRUE, skip=0, sep="\t"))
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
all(PED$Individual.ID == rownames(eigenvec)) == TRUE

# set colours
#BiocManager::install("RColorBrewer")
require("RColorBrewer")

# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
        "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
        "CLM","MXL","PEL","PUR",
        "CDX","CHB","CHS","JPT","KHV",
        "CEU","FIN","GBR","IBS","TSI",
        "BEB","GIH","ITU","PJL","STU"))

col <- colorRampPalette(c(
        "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
        "forestgreen","forestgreen","forestgreen","forestgreen",
        "grey","grey","grey","grey","grey",
        "royalblue","royalblue","royalblue","royalblue","royalblue",
        "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

# generate PCA bi-plots
project.pca <- eigenvec
summary(project.pca)


par(mar=c(5,5,5,5), cex=2.0, cex.main=7, cex.axis=2.75, cex.lab=2.75, mfrow=c(1,2))

plot(project.pca[,1], project.pca[,2], type="n", main="A", adj=0.5, xlab="First component", ylab="Second component", font=2, font.lab=2)
points(project.pca[,1], project.pca[,2], col=col, pch=20, cex=2.25)
legend("bottomright", bty="n", cex=3.0, title="", c("African","Hispanic","East-Asian","Caucasian","South Asian"), fill=c("yellow","forestgreen","grey","royalblue","black"))

plot(project.pca[,1], project.pca[,3], type="n", main="B", adj=0.5, xlab="First component", ylab="Third component", font=2, font.lab=2)
points(project.pca[,1], project.pca[,3], col=col, pch=20, cex=2.25)
