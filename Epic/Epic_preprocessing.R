setwd("~/Epic/raw_data/")
taitung_simon_bim = read.table("./taitung_simon_merge.bim", stringsAsFactors = FALSE)
str(taitung_simon_bim)
hgdp_bim = read.table("./HGDP_Map.bim", stringsAsFactors = FALSE)
str(hgdp_bim)
t_merge = merge(taitung_simon_bim, hgdp_bim, by.x = "V2", by.y = "V2")
str(t_merge)
write.table(t_merge[,1], file = "./merge_snp_list1.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile taitung_simon_merge --extract merge_snp_list1.txt --make-bed --out taitung_simon_extract")
system("plink --bfile HGDP_Map --extract merge_snp_list1.txt --make-bed --out hgdp_extract")
system("plink --bfile taitung_simon_extract --bmerge hgdp_extract --make-bed --out test1")
system("plink --bfile taitung_simon_extract --flip test1-merge.missnp --make-bed --out taitung_simon_extract_flip")
system("plink --bfile taitung_simon_extract_flip --bmerge hgdp_extract --make-bed --out test2")
system("plink --bfile taitung_simon_extract_flip --exclude test2-merge.missnp --make-bed --out exclude_taitung_simon_ex_flip")
system("plink --bfile exclude_taitung_simon_ex_flip --bmerge hgdp_extract --make-bed --out test3")
system("mv test3.bed taitung_simon_hgdp.bed")
system("mv test3.bim taitung_simon_hgdp.bim")
system("mv test3.fam taitung_simon_hgdp.fam")

sea_bim = read.table("./SEA_730K.bim", stringsAsFactors = FALSE)
str(sea_bim)
tai_simon_hgdp_bim = read.table("./taitung_simon_hgdp.bim", stringsAsFactors = FALSE)
str(tai_simon_hgdp_bim)
t_merge = merge(tai_simon_hgdp_bim, sea_bim, by.x = "V2", by.y = "V2")
str(t_merge)
write.table(t_merge[,1], file = "./merge_snp_list2.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile taitung_simon_hgdp --extract merge_snp_list2.txt --make-bed --out taitung_simon_hgdp_extract")
system("plink --bfile SEA_730K --extract merge_snp_list2.txt --make-bed --out sea_extract")
system("plink --bfile taitung_simon_hgdp_extract --bmerge sea_extract --make-bed --out test1")
system("plink --bfile taitung_simon_hgdp_extract --flip test1-merge.missnp --make-bed --out taitung_simon_hgdp_extract_flip")
system("plink --bfile taitung_simon_hgdp_extract_flip --bmerge sea_extract --make-bed --out test2")
system("mv test2.bed epic.bed")
system("mv test2.bim epic.bim")
system("mv test2.fam epic.fam")

epic_fam = read.table("./epic.fam", stringsAsFactors = FALSE)
str(epic_fam)
epic_fam$pop = NA

simon_sample_id = read.csv("~/Simon/raw_data/nature18964-s2.csv", header = TRUE, stringsAsFactors = FALSE)
str(simon_sample_id)
epic_fam$pop[!is.na(match(epic_fam[,2], simon_sample_id[,3]))] = simon_sample_id[,9][na.omit(match(epic_fam[,2], simon_sample_id[,3]))]

taitung_list = read.csv("~/PCA/R_analysis/96samples_sex_list.csv", header = F, stringsAsFactors = FALSE)
str(taitung_list)
epic_fam$pop[!is.na(match(epic_fam[,2], taitung_list[,2]))] = taitung_list[,4][na.omit(match(epic_fam[,2], taitung_list[,2]))]

hgdp_list = read.table("~/PCA/R_analysis/HGDP-country-list.txt", header = T, stringsAsFactors = F)
str(hgdp_list)
parse = unlist(strsplit(hgdp_list$ID_Pop,"_"))
hgdp_list["HGDP_ID"] = parse[seq(1,length(parse),2)]
hgdp_list["pop"] = parse[seq(2,length(parse),2)]
hgdp_list$pop = gsub("\\s","",hgdp_list$pop)
str(hgdp_list)
epic_fam$pop[!is.na(match(epic_fam[,2], hgdp_list[,2]))] = hgdp_list[,3][na.omit(match(epic_fam[,2], hgdp_list[,2]))]

sea_list = read.csv("~/SoutheastAsia/SEA_ind_info.csv", header = FALSE, stringsAsFactors = FALSE)
str(sea_list)
epic_fam$pop[!is.na(match(epic_fam[,2], sea_list[,1]))] = sea_list[,2][na.omit(match(epic_fam[,2], sea_list[,1]))]

write.table(epic_fam[,c(1,2,7)], file = "./data.clust", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(table(epic_fam$pop), file = "./count_table.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

epic_fam$pop[epic_fam$pop %in% "Ami"] = "Amis"
epic_fam$pop[epic_fam$pop %in% "Igorot"] = "Kankanaey"
epic_fam$pop[epic_fam$pop %in% "Kinh"] = "Vietnamese"
epic_fam$pop[epic_fam$pop %in% "Papuan"] = "Papuans"
pop_list = c("Amis", "Atayal", "Australian", "Bunun", "Burmese", "Dunsun", "Filipino", "Hawaiian",
             "Kankanaey", "Kinh", "Malay", "Maori", "Murut", "Paiwan", "Papuans", "Puyuma", "Vietnamese")
write.table(pop_list, file = "./asia_pop.list", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("epic", "~/Epic/raw_data/data.clust", "~/Epic/raw_data/asia_pop.list", "epic_selct")

epic_fam[epic_fam$V2 == "TDC480",]
epic_fam[epic_fam$V2 == "TDC480_2",]
epic_fam = epic_fam[-c(88,95),]
write.table(epic_fam[,c(1,2,7)], file = "./data.clust2", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("epic", "~/Epic/raw_data/data.clust2", "~/Epic/raw_data/asia_pop.list", "epic_select_no_outliner")
table(epic_fam$pop)

########################### Merge with Pan-Asia ################################
pan_asia_bim = read.table("./Genotypes_All.bim", stringsAsFactors = FALSE)
str(pan_asia_bim)
epic_bim = read.table("./epic.bim", stringsAsFactors = FALSE)
str(epic_bim)
t_merge = merge(pan_asia_bim, epic_bim, by.x = "V2", by.y = "V2")
str(t_merge)
write.table(t_merge[,1], file = "./epic_pan_merge_list1.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile epic --extract epic_pan_merge_list1.txt --make-bed --out epic_extract")
system("plink --bfile Genotypes_All --extract epic_pan_merge_list1.txt --make-bed --out pan_extract")
system("plink --bfile epic_extract --bmerge pan_extract --make-bed --out test3")
system("plink --bfile epic_extract --exclude test3-merge.missnp --make-bed --out exclude_epic")
system("plink --bfile exclude_epic --bmerge pan_extract --make-bed --out test4")
system("mv test4.bed epic_pan.bed")
system("mv test4.bim epic_pan.bim")
system("mv test4.fam epic_pan.fam")

epic_pan_fam = read.table("./epic_pan.fam", stringsAsFactors = FALSE)
str(epic_pan_fam)
epic_data_clust = read.table("./data.clust", stringsAsFactors = FALSE)
str(epic_data_clust)
epic_pan_fam$pop = NA

# Fill in Pan-Asia population list
library(stringr)
for(i in 1:nrow(epic_pan_fam)){
    pop =  str_sub(epic_pan_fam$V2[i],start = 1,end = 5)
    epic_pan_fam$pop[i] = pop 
}
epic_pan_fam$pop[!is.na(match(epic_pan_fam[,2], epic_data_clust[,2]))] = epic_data_clust[,3][na.omit(match(epic_pan_fam[,2], epic_data_clust[,2]))]
# Revice several population names
for(p in 1:nrow(epic_pan_fam)){
    if(epic_pan_fam$pop[p] == "CEU-N" | 
       epic_pan_fam$pop[p] == "CHB-N" |
       epic_pan_fam$pop[p] == "JPT-N" |
       epic_pan_fam$pop[p] == "YRI-N"){
        epic_pan_fam$pop[p] = substr(epic_pan_fam$pop[p],1,3)
    }
}
epic_pan_fam$pop[epic_pan_fam$pop %in% "Ami"] = "Amis"
epic_pan_fam$pop[epic_pan_fam$pop %in% "Igorot"] = "Kankanaey"
epic_pan_fam$pop[epic_pan_fam$pop %in% "Kinh"] = "Vietnamese"
epic_pan_fam$pop[epic_pan_fam$pop %in% "Papuan"] = "Papuans"
epic_pan_fam$pop[epic_pan_fam$pop %in% "AX-AT"] = "Atayal"
epic_pan_fam$pop[epic_pan_fam$pop %in% "AX-AM"] = "Amis"
epic_pan_fam[epic_pan_fam$V2 == "TDC480",]
epic_pan_fam[epic_pan_fam$V2 == "TDC480_2",]
epic_pan_fam = epic_pan_fam[-c(88,95),]
write.table(epic_pan_fam[,c(1,2,7)], file = "./data.clust3", quote = FALSE, row.names = FALSE, col.names = FALSE)
pop_list = c("Amis", "Atayal", "Australian","Bunun", "Burmese", "Dunsun", "Filipino", "Hawaiian",
             "Kankanaey", "Kinh", "Malay", "Murut", "Paiwan", "Papuans", "Puyuma",
             "ID-ML", "ID-DY", "ID-TR", "PI-IR", "PI-UB", "PI-MT", "PI-UN", "PI-UI", "PI-AG",
             "PI-AE", "PI-AT", "PI-MW")
write.table(pop_list, file = "./pan_pop.list", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("epic_pan", "~/Epic/raw_data/data.clust3", "~/Epic/raw_data/pan_pop.list", "epic_pan_select")

# Remove Australian, Papuans
pop_list = c("Amis", "Atayal","Bunun", "Burmese", "Dunsun", "Filipino", "Hawaiian",
             "Kankanaey", "Kinh", "Malay", "Murut", "Paiwan", "Puyuma",
             "ID-ML", "ID-DY", "ID-TR", "PI-IR", "PI-UB", "PI-MT", "PI-UN", "PI-UI", "PI-AG",
             "PI-AE", "PI-AT", "PI-MW")
write.table(pop_list, file = "./pan_pop.list2", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("epic_pan", "~/Epic/raw_data/data.clust3", "~/Epic/raw_data/pan_pop.list2", "epic_pan_select2")

# Remove PI-MW, PI-AT, PI-AE, PI-AG, PI-IR
pop_list = c("Amis", "Atayal","Bunun", "Burmese", "Dunsun", "Filipino", "Hawaiian",
             "Kankanaey", "Kinh", "Malay", "Murut", "Paiwan", "Puyuma",
             "ID-ML", "ID-DY", "ID-TR", "PI-UB", "PI-MT", "PI-UN", "PI-UI")
write.table(pop_list, file = "./pan_pop.list3", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("epic_pan", "~/Epic/raw_data/data.clust3", "~/Epic/raw_data/pan_pop.list3", "epic_pan_select3")

# Remove PI-MW, PI-AT, PI-AE, PI-AG, PI-IR ADD PI-MA
pop_list = c("Amis", "Atayal","Bunun", "Burmese", "Dunsun", "Filipino", "Hawaiian",
             "Kankanaey", "Kinh", "Malay", "Murut", "Paiwan", "Puyuma",
             "ID-ML", "ID-DY", "ID-TR", "PI-UB", "PI-MT", "PI-UN", "PI-UI", "PI-MA")
write.table(pop_list, file = "./pan_pop.list4", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("epic_pan", "~/Epic/raw_data/data.clust3", "~/Epic/raw_data/pan_pop.list4", "epic_pan_select4")

################################# TREEMIX SECTION ##############################
setwd("~/Epic/TreeMix/")
system("plink --bfile epic_pan --freq --missing --within data.clust3 --out epic_pan")
system("gzip -c epic_pan.frq.strat > epic_pan.frq.strat.gz")
system("python plink2treemix.py epic_pan.frq.strat.gz epic_pan.treemix.gz")

## IN SERVER
# treemix -i epic_pan.treemix.gz -root San -o epic_pan

n = "epic_pan"
file.temp = paste0("treemix_bootstrap_",n, ".sh")
line.0 = paste0("cd /home/chunyo/Epic")
line.1 = paste0("treemix -i epic_pan.treemix.gz -root San -bootstrap -k 15 -o epic_pan")
cat(line.0, file= file.temp , append=TRUE,sep="\n")
cat(line.1, file= file.temp , append=TRUE,sep="\n")

file.temp = paste0("pbs_treemix_bootstrap_",n,".sh")
line.1 = paste0("#PBS -N treemix_bootstrap_",n)
line.2 = paste0("#PBS -e treemix_bootstrap_",n,".err")
line.3 = paste0("#PBS -o treemix_bootstrap_",n,".out")
line.4 = paste0("cd $PBS_O_WORKDIR")
line.5 = paste0("sh treemix_bootstrap_",n,".sh")
cat(line.1, file= file.temp , append=TRUE,sep="\n")
cat(line.2, file= file.temp , append=TRUE,sep="\n")
cat(line.3, file= file.temp , append=TRUE,sep="\n")
cat(line.4, file= file.temp , append=TRUE,sep="\n")
cat(line.5, file= file.temp , append=TRUE,sep="\n")

# qsub pbs_treemix_bootstrap_epic_pan.sh
source("~/TreeMix/plotting_funcs.R")
plot_tree("epic_pan")
