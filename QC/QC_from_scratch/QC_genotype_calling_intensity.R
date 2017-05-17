### Genotype Calling and Signal Intensities ###

setwd("/Users/chenchun-yu/Rstudio/QC_GWAS")
source("http://bioconductor.org/biocLite.R")
#Install package "RSQLite" to liberary packages_LB.
biocLite()
install.packages("RSQLite",dependencies = T,lib="./packages_LB/")
#Loading the library.
library(RSQLite)


### Find the file that contains "X,Y values". ###
#Creat a empty database called "sample1DB".
con=dbConnect(dbDriver("SQLite"),dbname="sample1DB")
#Populate the database by writting the table from file "8samples_FinalReport1.txt".
dbWriteTable(con,"sample1_snps","8samples_FinalReport1.txt",header=T,append=T,skip=10,sep="\t")
#Select the columns 1,3,4,18,19.
##Column 1 contains "SNP.Name".
##Column 3 and 4 contain "Allele1,2...Top".
##Column 18 and 19 contain "X,Y values" which represent intensity.
snp=dbGetQuery(con,"select * from sample1_snps")[,c(1,3,4,18,19)]


### Find the file that contains "allele genotype" in A/B format. ### 
#Creat a empty SQL database called "SNP_Map".
snp_mapDB=dbConnect(dbDriver("SQLite"),dbname="SNP_Map")
#Populate the database by writting the table from file "SNP Map.txt".
dbWriteTable(snp_mapDB,"snp_map","SNP Map.txt",append=T,header=T,sep="\t")
#Column 11 contains the genotype of sample "TDC480". 
TDC480GType=dbGetQuery(snp_mapDB,"select * from snp_map")[11]
TDC480GType=as.factor(TDC480GType$TDC480.GType)
#Add sample "TDC480" genotype to "snp" table.
snp$genotype <- TDC480GType


### Make a X/Y intensity plot for all samples. ###
plot(snp$X,snp$Y,col=snp$genotype,pch=as.numeric(snp$genotype),
     xlab = "x",ylab = "y",main = "sample_one",cex.main=0.9)
legend("topright",paste(levels(snp$genotype),"(",summary(snp$genotype),")",sep = ""),
       col = 1:length(levels(snp$genotype)),pch=1:length(levels(snp$genotype)),cex=0.7)


### Make a X/Y intensity plot for firt 2500 samples. ###
plot(head(snp$X,n=2500),head(snp$Y,n=2500),col=head(snp$genotype,n=2500),pch=as.numeric(head(snp$genotype,n=2500)),
     xlab = "x",ylab = "y",main = "sample_one",cex.main=0.9)
legend("topright",paste(levels(head(snp$genotype,n=2500)),
      "(",summary(head(snp$genotype,n=2500)),")",sep = ""),
       col = 1:length(levels(head(snp$genotype,n=2500))),
       pch=1:length(levels(head(snp$genotype,n=2500))),cex=0.7)