setwd("/Users/chenchun-yu/Rstudio/QC_GWAS")
source("http://bioconductor.org/biocLite.R")
#Install package "RSQLite" to liberary packages_LB.
biocLite()
install.packages("RSQLite",dependencies = T,lib="./packages_LB/")
#Loading the library.
library(RSQLite)

#Creat a empty SQL database called "SNP_Map".
snp_mapDB=dbConnect(dbDriver("SQLite"),dbname="SNP_Map")
#Populate the database by writting the table from file "SNP Map.txt".
dbWriteTable(snp_mapDB,"snp_map","SNP Map.txt",append=T,header=T,sep="\t")

#List table of the database.
dbListTables(snp_mapDB)
#List fields of the table "snp_map" in database "snp_mapDB".
dbListFields(snp_mapDB,"snp_map")
#List first three rows of the table.
head(dbReadTable(snp_mapDB,"snp_map"),n=3)
#Get the nunber of total SNP in this case.
dbGetQuery(snp_mapDB,"select count (*) from snp_map")

#Select the columns which contain genotype.
##Column 11 contains the genotype of sample "TDC480". 
TDC480GType=dbGetQuery(snp_mapDB,"select * from snp_map")[11]
##Column 46 contains the genotype of sample "TDC480_2".
TDC480.2GType=dbGetQuery(snp_mapDB,"select * from snp_map")[46]

#Find the different genotype between TDC480 and TDC480_2
different=which((TDC480GType != TDC480.2GType))
#Make a table called "compare_table" to record the differences.
compare_table=data.frame(TDC480=TDC480GType[different,],TDC480_2=TDC480.2GType[different,])

#Added the columns to "compare_table" which combined two genotypes.
compare_table$type <- as.factor(paste(compare_table$TDC480,compare_table$TDC480_2,sep = "/"))
#Show the total number for each unique genotypes.
summary(compare_table$type)
#Make a barplot.
barplot(summary(compare_table$type))