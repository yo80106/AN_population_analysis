setwd("/Users/chenchun-yu/Rstudio/QC_GWAS")
source("http://bioconductor.org/biocLite.R")
#Install package "RSQLite" to liberary packages_LB.
biocLite()
install.packages("RSQLite",dependencies = T,lib="./packages_LB/")
#Loading the library.
library(RSQLite)

#Creat a empty database called "sample1DB".
con=dbConnect(dbDriver("SQLite"),dbname="sample1DB")
#Populate the database by writting the table from file "8samples_FinalReport1.txt".
dbWriteTable(con,"sample1_snps","8samples_FinalReport1.txt",header=T,append=T,skip=10,sep="\t")

#List table of the database. 
dbListTables(con)
#List fields of the table "sample1_snps" in database "con".
dbListFields(con,"sample1_snps")
#List first three rows of the table.
head(dbReadTable(con, "sample1_snps"),n=3)

#Get the nunber of total SNP in this case.
dbGetQuery(con,"select count (*) from sample1_snps")
#Get the numbers of sample in this table.
sample_count=dbGetQuery(con,"select distinct 'Sample.ID' from sample1_snps")
dim(sample_count)

#Select the columns 1,3,4,18,19.
##Column 1 contains "SNP.Name".
##Column 3 and 4 contain "Allele1,2...Top".
##Column 18 and 19 contain "X,Y values" which represent intensity.
snp=dbGetQuery(con,"select * from sample1_snps")[,c(1,3,4,18,19)]
dim(snp)
#Get the numbers of missing values(-), ATCG, and indels(D/T).
snp$Allele1...Top=factor(snp$Allele1...Top)
snp$Allele2...Top=factor(snp$Allele2...Top)
summary(snp$Allele1...Top)
summary(snp$Allele2...Top)
#Find out the unique combinition of Allele1 and Allele2 genotype.
unique(paste(snp$Allele1...Top,snp$Allele2...Top))