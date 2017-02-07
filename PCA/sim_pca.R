source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")
install.packages("scatterplot3d")
library(gdsfmt)
library(SNPRelate)
library(scatterplot3d)
genofile <- snpgdsOpen(snpgdsExampleFileName())
(snp <- read.gdsn(index.gdsn(genofile, "snp.allele")))
(geno <- read.gdsn(index.gdsn(genofile, "genotype")))
(info <- read.gdsn(index.gdsn(genofile, "sample.annot")))
str(snp)
str(geno)
class(geno)
dim(geno)

pca <- snpgdsPCA(genofile)
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

geno_cov = cov(geno[,1:9088])
heptathlon_pca <- prcomp(geno, scale = TRUE)
summary(geno_pcacov, loadings = TRUE)



sim = geno[1:12,1:10]
sim[1,]
sim[,1]
colMeans(sim)
sim_cov = cov(sim[,1:100])
# sim_cor = cor(sim[,1:100])
scatterplot3d(sim[,1],sim[,2],sim[,3])


pcacov = princomp(covmat=sim_cov)
sim_pca = prcomp(sim, scale = TRUE)
summary(pcacov, loadings = TRUE)

a1 = pcacov$scale
drop(scale(hm, center = center, 
           scale = scale) %*% 
         heptathlon_pca$rotation[,1])

plot(pcacor$sdev^2, 
     xlab = "Component number",
     ylab = "Component variance", 
     type = "l", main = "Scree diagram")