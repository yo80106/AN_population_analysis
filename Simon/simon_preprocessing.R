setwd("~/Simon/")
simon_sample_id = read.csv("./raw_data/nature18964-s2.csv", header = TRUE, stringsAsFactors = FALSE)
str(simon_sample_id)

oceania_ind = simon_sample_id$Sample.ID..Illumina.[simon_sample_id$Region == "Oceania"]
# I found that LP6005443-DNA_G07 has some problem, so I remove it manually
oceania_ind = oceania_ind[-20]
simon_sample_id[simon_sample_id$Region == "Oceania",][,c(3,9,10)]
system("mkdir Oceania")
setwd("~/Simon/raw_data/")
for(i in oceania_ind){
    system(paste0("plink -vcf ",i,".annotated.nh2.variants.vcf.gz --exclude dot_rm.txt --double-id --make-bed --out ./Oceania/",i))
}
setwd("~/Simon/raw_data/Oceania/")
merge_bfile = list.files("./",pattern = "(.bim|.fam|.bed)$")
t_matrix = matrix(merge_bfile, ncol = 3, byrow = TRUE)
write.table(t_matrix[-1,],"./merge_list.txt",col.names = F,row.names = F,quote = F)
t_matrix[1,]
system("plink --bfile LP6005441-DNA_A03 --merge-list merge_list.txt --make-bed --out merge1")
for(z in oceania_ind){
    system(paste0("plink --bfile ",z," --exclude merge1-merge.missnp --make-bed --out ",z))
}
system("plink --bfile LP6005441-DNA_A03 --merge-list merge_list.txt --make-bed --out merge2")

######################### ALL DATA IN SIMON PROJECT ############################
system("mkdir bfile")
for(i in simon_sample_id[,3]){
    system(paste0("plink -vcf ",i,".annotated.nh2.variants.vcf.gz --exclude dot_rm.txt --double-id --make-bed --out ./bfile/",i))
}
a = list.files("./", pattern = ".gz$")
# A dot (\\.), then any character (.) any number of times (*) until the end of the string ($).
a = gsub("\\..*$", "", a)
b = simon_sample_id$Sample.ID..Illumina.[simon_sample_id$Embargo.level..X.Fully.Public..Y.Signed.Letter. == "Y"]
c = setdiff(simon_sample_id[,3], a)
setdiff(c,b)

setwd("~/Simon/raw_data/bfile/")
simon_bfile = list.files("./",pattern = "(.bim|.fam|.bed)$")
t_matrix = matrix(simon_bfile, ncol = 3, byrow = TRUE)
write.table(t_matrix[-1,],"./merge_list.txt",col.names = F,row.names = F,quote = F)
# I found that LP6005443-DNA_G07 has some problem, so I remove it manually
t_matrix[1,]
system("plink --bfile LP6005441-DNA_A01 --merge-list merge_list.txt --make-bed --out merge1")
new_merge_list = read.table("./merge_list.txt")
new_id = gsub("\\..*$", "",new_merge_list[,1])
for(z in new_id){
    system(paste0("plink --bfile ",z," --exclude merge1-merge.missnp --make-bed --out ",z))
}
system("plink --bfile LP6005441-DNA_A01 --merge-list merge_list.txt --make-bed --out merge2")
## IN SERVER
simon_bim = read.table("./merge2.bim", stringsAsFactors = FALSE)
str(simon_bim)
taitung_bim = read.table("./new_clean_per_snp.bim", header = FALSE, stringsAsFactors = FALSE)
str(taitung_bim)
t_merge = merge(taitung_bim, simon_bim, by.x = "V2", by.y = "V2")
dim(t_merge)
head(t_merge)
write.table(t_merge[,1], file = "./merge_snp_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile new_clean_per_snp --extract merge_snp_list.txt --make-bed --out taitung_extract")
system("plink --bfile merge2 --extract merge_snp_list.txt --make-bed --out simon_extract")
system("plink --bfile taitung_extract --bmerge simon_extract --make-bed --out test1")
system("plink --bfile taitung_extract --flip test1-merge.missnp --make-bed --out taitung_extract_flip")
system("plink --bfile taitung_extract_flip --bmerge simon_extract --make-bed --out test2")
system("plink --bfile taitung_extract_flip --exclude test2-merge.missnp --make-bed --out exclude_taitung_ex_flip")
system("plink --bfile exclude_taitung_ex_flip --bmerge simon_extract --make-bed --out test3")
system("mv test3.bed taitung_simon_merge.bed")
system("mv test3.bim taitung_simon_merge.bim")
system("mv test3.fam taitung_simon_merge.fam")

merge_fam = read.table("./taitung_simon_merge.fam", stringsAsFactors = FALSE)
source("~/Git/PCA/pca_function.R")
taitung_list = read.csv("~/PCA/R_analysis/96samples_sex_list.csv", header = F, stringsAsFactors = FALSE)
merge_fam$pop = NA
merge_fam$pop[!is.na(match(merge_fam[,2], taitung_list[,2]))] = taitung_list[,4][na.omit(match(merge_fam[,2], taitung_list[,2]))]
merge_fam$pop[!is.na(match(merge_fam[,2], simon_sample_id[,3]))] = simon_sample_id[,9][na.omit(match(merge_fam[,2], simon_sample_id[,3]))]

data.clust1 = write.table(merge_fam[,c(1,2,7)], file = "./data.clust1", quote = FALSE, col.names = FALSE, row.names = FALSE)
pop.list1 = write.table(unique(merge_fam$pop), file = "./pop.list1", quote = FALSE, col.names = FALSE, row.names = FALSE)
my_pca("taitung_simon_merge", "~/Simon/raw_data/bfile/data.clust1", "~/Simon/raw_data/bfile/pop.list1", "taitung_simon_merge")

## Asia, Oceania PCA
asia_draft_id = simon_sample_id[simon_sample_id$Region %in% c("EastAsia", "Oceania", "SouthAsia"),][,c(3,7,9,10)]
asia_id = asia_draft_id[asia_draft_id$Embargo.level..X.Fully.Public..Y.Signed.Letter. %in% "X",]
write.table(unique(asia_id$Population.ID), file = "./asia_id.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
my_pca("taitung_simon_merge", "~/Simon/raw_data/bfile/data.clust1", "~/Simon/raw_data/bfile/asia_id.txt", "taitung_simon_asia")