setwd("~/PCA/R_analysis/HGDP_Pan_merge/NO_LD")
############# Change RS Number ##############
system("plink --bfile new_clean_per_snp --update-map ~/PCA/R_analysis/rsids.txt --update-name --make-bed --out new_clean_per_snp_changed")
## Check by R
rs_list = read.table("~/PCA/R_analysis/rsids.txt",header = TRUE,na.strings = ".")
taitung_bim = read.table("~/PCA/R_analysis/HGDP_Pan_merge/NO_LD/new_clean_per_snp.bim")
merge_id = merge(taitung_bim,rs_list,by.x = "V2",by.y = "Name")
i = sapply(merge_id, is.factor)
merge_id[i] = lapply(merge_id[i], as.character)
na = sum(is.na(merge_id$RsID))
nrow(merge_id)-na
############# HGDP and Pan-Asia Merge ##############
hgdp_bim = read.table("~/PCA/R_analysis/HGDP_Pan_merge/HGDP_Map.bim",stringsAsFactors = FALSE)
pan_bim = read.table("~/PCA/R_analysis/HGDP_Pan_merge/Genotypes_All.bim",stringsAsFactors = FALSE)
snplist = merge(hgdp_bim,pan_bim,by = "V2")
write.table(snplist$V2,file = "~/PCA/R_analysis/HGDP_Pan_merge/combine_snp_list.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
############# PLINK SECTION ##############
system("plink --bfile Genotypes_All --extract combine_snp_list.txt --make-bed --out extract_pan")
system("plink --bfile HGDP_Map --extract combine_snp_list.txt --make-bed --out extract_hgdp")
system("plink --bfile extract_pan --bmerge extract_hgdp.bed extract_hgdp.bim extract_hgdp.fam --make-bed --out combine")
system("mv combine-merge.missnp flipsnp.txt")
system("plink --bfile extract_pan --flip flipsnp.txt --make-bed --out flip-extract_pan")
system("plink --bfile flip-extract_pan --bmerge extract_hgdp.bed extract_hgdp.bim extract_hgdp.fam --make-bed --out flip_combine") 
############# HGDP, Pan-Asia and Taitung Merge ##############
flip_combine_bim = read.table("~/PCA/R_analysis/HGDP_Pan_merge/flip_combine.bim",stringsAsFactors = FALSE)
taitung_bim = read.table("~/PCA/R_analysis/HGDP_Pan_merge/NO_LD/new_clean_per_snp_changed.bim",stringsAsFactors = FALSE) 
snplist2 = merge(flip_combine_bim,taitung_bim,by = "V2")
write.table(unique(snplist2$V2),file = "~/PCA/R_analysis/HGDP_Pan_merge/NO_LD/combine_taitung_snp_list.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
############# PLINK SECTION ##############
system("plink --bfile ../flip_combine --extract combine_taitung_snp_list.txt --make-bed --out flip_combine_ex_taitung")
system("plink --bfile new_clean_per_snp_changed --extract combine_taitung_snp_list.txt --make-bed --out extract_taitung")
system("plink --bfile flip_combine_ex_taitung --bmerge extract_taitung.bed extract_taitung.bim extract_taitung.fam --make-bed --out taitung_combine_noflip")
system("mv taitung_combine_noflip-merge.missnp flip-taitung-combine-list.txt")
system("plink --bfile flip_combine_ex_taitung --flip flip-taitung-combine-list.txt --make-bed --out flip-flip_combine")
system("plink --bfile flip-flip_combine --bmerge extract_taitung.bed extract_taitung.bim extract_taitung.fam --make-bed --out flip_taitung_combine")

## Venn Diagram
install.packages('VennDiagram')
library(VennDiagram)
nhgdp_bim = nrow(hgdp_bim)
npan_bim = nrow(pan_bim)
ntaitung_bim = nrow(taitung_bim)
a = merge(hgdp_bim,taitung_bim,by = "V2")
n12 = length(unique(a$V2))
b = merge(hgdp_bim,pan_bim,by = "V2")
n23 = length(unique(b$V2))
c = merge(pan_bim,taitung_bim,by = "V2")
n13 = length(unique(c$V2))
n123 = length(unique(snplist2$V2))
pdf("~/PCA/R_analysis/venn.pdf")
draw.triple.venn(area1 = ntaitung_bim, area2 = nhgdp_bim, area3 = npan_bim, 
                 n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = c("Taitung","HGDP","Pan-Asia"),
                 lty = "blank",fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()