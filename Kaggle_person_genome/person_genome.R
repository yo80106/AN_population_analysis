setwd("~/Kaggle/genome/")
raw_data = read.csv("./genome_zeeshan_usmani.csv", header = T, stringsAsFactors = F)
write.table(raw_data, file = "./zeeshan_usmani.txt", quote = F, row.names = F, col.names = F, sep = "\t")
system("plink --23file zeeshan_usmani.txt Zeeshan Usmani --out zeeshan_usmani")

zee_bim = read.table("./zeeshan_usmani.bim", stringsAsFactors = FALSE)
str(zee_bim)
hgdp_bim = read.table("./HGDP_Map.bim", stringsAsFactors = FALSE)
str(hgdp_bim)
t_merge = merge(zee_bim, hgdp_bim, by.x = "V2", by.y = "V2")
str(t_merge)
write.table(t_merge[,1], file = "./merge_snp_list1.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile zeeshan_usmani --extract merge_snp_list1.txt --make-bed --out zeeshan_usmani_extract")
system("plink --bfile HGDP_Map --extract merge_snp_list1.txt --make-bed --out hgdp_extract")
system("plink --bfile zeeshan_usmani_extract --bmerge hgdp_extract --make-bed --out test1")
system("plink --bfile zeeshan_usmani_extract --flip test1-merge.missnp --make-bed --out zeeshan_usmani_extract_flip")
system("plink --bfile zeeshan_usmani_extract_flip --bmerge hgdp_extract --make-bed --out test2")
system("mv test2.bed zeeshan_usmani_hgdp.bed")
system("mv test2.bim zeeshan_usmani_hgdp.bim")
system("mv test2.fam zeeshan_usmani_hgdp.fam")

hgdp_list = read.table("~/PCA/R_analysis/HGDP-country-list.txt", header = T, stringsAsFactors = F)
str(hgdp_list)
parse = unlist(strsplit(hgdp_list$ID_Pop,"_"))
hgdp_list["HGDP_ID"] = parse[seq(1,length(parse),2)]
hgdp_list["pop"] = parse[seq(2,length(parse),2)]
hgdp_list$pop = gsub("\\s","",hgdp_list$pop)
str(hgdp_list)

zee_fam = read.table("./zeeshan_usmani_hgdp.fam", stringsAsFactors = FALSE)
str(zee_fam)
zee_fam$pop = NA
zee_fam$pop[!is.na(match(zee_fam[,2], hgdp_list[,2]))] = hgdp_list[,3][na.omit(match(zee_fam[,2], hgdp_list[,2]))]
zee_fam$pop[1044] = "Zeeshan"
write.table(zee_fam[,c(1,2,7)], file = "./data.clust1", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(unique(zee_fam$pop), file = "./pop.list1", quote = F, row.names = F, col.names = F)
my_pca("zeeshan_usmani_hgdp", "~/Kaggle/genome/data.clust1", "~/Kaggle/genome/pop.list1", "zee_hgdp_1")

pop_list = c("Brahui","Balochi", "Hazara", "Makrani", "Sindhi", "Pathan", "Kalash", "Burusho", "Druze", "Bedouin", "Zeeshan")
write.table(pop_list, file = "./pop.list2", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("zeeshan_usmani_hgdp", "~/Kaggle/genome/data.clust1", "~/Kaggle/genome/pop.list2", "zee_hgdp_2")

pop_list = c("Brahui","Balochi", "Hazara", "Makrani", "Sindhi", "Pathan", "Burusho", "Zeeshan")
write.table(pop_list, file = "./pop.list3", quote = FALSE, row.names = FALSE, col.names = FALSE)
my_pca("zeeshan_usmani_hgdp", "~/Kaggle/genome/data.clust1", "~/Kaggle/genome/pop.list3", "zee_hgdp_3")