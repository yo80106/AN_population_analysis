setwd("~")
# Install packages for PCA analysis
source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")

library(gdsfmt)
library(SNPRelate)

# Exclude sample 18
fam = read.table("~/PCA/R_analysis/new_clean_per_snp.fam")
bim = read.table("~/PCA/R_analysis/new_clean_per_snp.bim")
ind_info = read.csv("~/PCA/R_analysis/96samples_sex_list.csv",header = F)
ind_info = ind_info[-18,]
ind_info = data.frame(ind_info,row.names = 1:95)
ind_info = droplevels(ind_info)
table(ind_info$V4)
# Check the individual ID between ind_info and fam
sum(ind_info$V2!=fam$V2)

# Conversion from PLINK BED to GDS
bedfile = "~/PCA/R_analysis/new_clean_per_snp.bed"
famfile = "~/PCA/R_analysis/new_clean_per_snp.fam"
bimfile = "~/PCA/R_analysis/new_clean_per_snp.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "taitung.gds")

genofile = snpgdsOpen("taitung.gds")
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))


pca = snpgdsPCA(genofile)
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
pop = ind_info$V4
tab = data.frame(sample.id = pca$sample.id,
                  pop = factor(pop)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
pdf("./taitung_pca.pdf")
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="PC 2", ylab="PC 1", main="PCA using all SNPs")
legend(-0.12,-0.35, legend=levels(tab$pop), pch="o", col=1:(nlevels(tab$pop)))
dev.off()
pdf("./taitung.pdf")
pc.percent = pca$varprop*100
lbls = paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)
dev.off()