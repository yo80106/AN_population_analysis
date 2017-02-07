########################## DATA PREPROCESSING ##################################
setwd("~/PCA/R_analysis/HGDP_Pan_merge/")
library(stringr)
all_fam = read.table("~/PCA/R_analysis/HGDP_Pan_merge/flip_taitung_combine.fam")
all_fam$pop = NA

# Fill in Pan-Asia population list
for(i in 1:nrow(all_fam)){
    pop =  str_sub(all_fam$V2[i],start = 1,end = 5)
    all_fam$pop[i] = pop 
}

# Fill in HGDP population list
hgdp_list = read.table("~/PCA/R_analysis/HGDP-country-list.txt",
                       header = T,stringsAsFactors = F)
parse = unlist(strsplit(hgdp_list$ID_Pop,"_"))
hgdp_list["HGDP_ID"] = parse[seq(1,length(parse),2)]
hgdp_list["pop"] = parse[seq(2,length(parse),2)]
hgdp_list$pop = gsub("\\s","",hgdp_list$pop)

for(t in 1:nrow(all_fam)){
    for(k in 1:nrow(hgdp_list)){
        if(all_fam$V2[t] == hgdp_list$HGDP_ID[k]){
            all_fam$pop[t] = hgdp_list$pop[k]
        }
    }
}

# Fill in Taitung population list
taitung_list = read.csv("~/PCA/R_analysis/96samples_sex_list.csv",header = F,
                        stringsAsFactors = FALSE)

for(t in 1:nrow(all_fam)){
    for(k in 1:nrow(taitung_list)){
        if(all_fam$V2[t] == taitung_list$V2[k]){
           all_fam$pop[t] = taitung_list$V4[k]
        }
    }
}

for(p in 1:nrow(all_fam)){
    if(all_fam$pop[p] == "CEU-N" | 
       all_fam$pop[p] == "CHB-N" |
       all_fam$pop[p] == "JPT-N" |
       all_fam$pop[p] == "YRI-N"){
        all_fam$pop[p] = substr(all_fam$pop[p],1,3)
    }
}

all_fam = all_fam[!(all_fam$pop %in% c("Han","Han−NChina","Miao","Japanese","CEU","JPT","YRI")),]
write.table(unique(all_fam$pop), file = "~/TreeMix/pop_list.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(all_fam[,c(1,2)], file = "~/PCA/R_analysis/HGDP_Pan_merge/final_merge_list.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)
system("plink --bfile flip_taitung_combine --keep final_merge_list.txt --make-bed --out final_merge")
# Remove high LD SNPs
system("plink --bfile final_merge --indep-pairwise 150 5 0.8 --out final_merge_no_LD")
system("plink --bfile final_merge --extract final_merge_no_LD.prune.in --make-bed --out final_merge_no_LD")
########################## Austronesian PCA SECTION ############################
pan_asia_language = read.delim("~/PCA/R_analysis/pan-asia-language.txt",header = TRUE,stringsAsFactors = FALSE)
austronesian = subset(pan_asia_language, Language.family == "  Austronesian")
austronesian_pop = gsub("\\s", "", austronesian$Pop.name)
austronesian_pop = c(austronesian_pop,"Paiwan","Amis","Bunun","Puyuma","TW-HA","TW-HB")
write.table(austronesian_pop, file = "~/Project/raw_data/austro_pop_list.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)

my_pca("final_merge_no_LD", "data.clust", "austro_pop_list.txt", "austronesian")
########################## Tree PCA SECTION ####################################
tree_pop = c("Paiwan","Amis","Bunun","Puyuma",
            "AX-AM","AX-AT","PI-UB","PI-UI",
            "PI-UN","ID-MT")
# tree_fam[tree_fam$V2 == "TDC480",]
# tree_fam[tree_fam$V2 == "TDC480_2",]
# tree_fam = tree_fam[-c(88,95),]
write.table(tree_pop,file = "~/Project/raw_data/tree_pop_list.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)
my_pca("final_merge_no_LD", "data.clust", "tree_pop_list.txt", "tree_pop")
########################## Philipines POPULATION PCA SECTION ###################
philipine_pop = c("Paiwan","Amis","Bunun","Puyuma",
                "AX-AM","AX-AT","PI-AG","PI-IR","PI-MA",
                "PI-MW","PI-UB","PI-UI","PI-UN")
write.table(philipine_pop,file = "~/Project/raw_data/philipine_pop_list.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)
my_pca("final_merge_no_LD", "data.clust", "philipine_pop_list.txt", "philipine_pop")
########################## TAIWAN POPULATION PCA SECTION #######################
pop_list = data.frame(c("Paiwan","Amis","Bunun","Puyuma",
                        "AX-AM","AX-AT","TW-HA","TW-HB"))
write.table(pop_list, file = "~/Project/raw_data/taiwan.txt", quote = F, row.names = F, col.names = F)
my_pca(plink.file = "final_merge_no_LD", data.clust = "data.clust", pop.list = "taiwan.txt", output.name = "taiwan_pop")
########################## FORMOSAN POPULATION PCA SECTION #####################
pop_list = data.frame(c("Paiwan","Amis","Bunun","Puyuma",
                        "AX-AM","AX-AT"))
write.table(pop_list, file = "~/Project/raw_data/aborigines.txt", quote = F, row.names = F, col.names = F)
my_pca(plink.file = "final_merge_no_LD", data.clust = "data.clust", pop.list = "aborigines.txt", output.name = "aborigines_pop")
########################## Han PCA SECTION #####################################
han_pop = c("CHB","TW-HB","TW-HA","SG-CN","CN-CC","CN-GA","CN-HM","CN-JI","CN-SH")
write.table(han_pop,file = "~/Project/raw_data/han_pop_list.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)
my_pca(plink.file = "final_merge_no_LD", data.clust = "data.clust", pop.list = "han_pop_list.txt", output.name = "han_pop")
system("plink --bfile final_merge --keep han_pop_list.txt --make-bed --out han")
###################### Populations admixed with Taiwanese aborigines ###########
pop_list = data.frame(c("Paiwan","Amis","Bunun","Puyuma",
                        "AX-AM","AX-AT","PI-MA",
                        "PI-UB","PI-UI","PI-UN",
                        "CN-CC","CN-JI","ID-MT"))
write.table(pop_list, file = "~/Project/raw_data/admixed.txt", quote = F, row.names = F, col.names = F)
my_pca(plink.file = "final_merge_no_LD", data.clust = "data.clust", pop.list = "admixed.txt", output.name = "admixed_with_pop")
########################## ALL MERGING DATA PCA SECTION ##############################
library(gdsfmt)
library(SNPRelate)

fam = read.table("~/PCA/R_analysis/HGDP_Pan_merge/final_merge.fam")
ind_info = all_fam
# Check the individual ID between ind_info and fam
sum(as.character(ind_info$V2)!=as.character(fam$V2))

# Conversion from PLINK BED to GDS
bedfile = "~/PCA/R_analysis/HGDP_Pan_merge/final_merge.bed"
famfile = "~/PCA/R_analysis/HGDP_Pan_merge/final_merge.fam"
bimfile = "~/PCA/R_analysis/HGDP_Pan_merge/final_merge.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "All.gds")

genofile = snpgdsOpen("All.gds")
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

pca = snpgdsPCA(genofile)
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
pop = ind_info$pop
tab = data.frame(sample.id = pca$sample.id,
                 pop = factor(pop)[match(pca$sample.id, sample.id)],
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)
# library(ggplot2)
# ggplot(tab, aes(EV2,EV1)) +
#     geom_point(aes(EV2,EV1,color = pop)) +
#     geom_text(aes(label=ifelse(pop == c("French","Yoruba"),as.character(pop),"")),hjust=0, vjust=0,color="black")
#     scale_colour_manual(breaks = df$zone, 
#                         values = unique(as.character(df$color.codes)))
      
#     theme(legend.direction ="vertical",legend.position = "right") +
#    guides(fill=guide_legend(ncol=2)) + 
#     scale_fill_discrete(lebels = levels(tab$pop))
# pdf("./global.pdf")
plot(tab$EV2, tab$EV1,col = rgb(214,214,214, maxColorValue = 255), xlab="PC 2", ylab="PC 1", main="PCA using all SNPs")
pca_pop = c("Paiwan","Amis","Bunun","Puyuma","AX-AM","AX-AT","AX-ME","Yoruba","Maya","French",
        "Palestinian","IN-WL","Sindhi","Papuan","PI-AE","PI-AG","ID-AL","ID-DY","TW-HA","TW-HB",
        "CHB","JP-ML","JP-RK","CN-HM","TH-HM")
id_nrow = which(tab$pop %in% pca_pop)
id_row_pop = tab$pop[tab$pop %in% pca_pop]
colors = c(rgb(170,93,152, maxColorValue=255),
      rgb(1,157,249, maxColorValue=255),
      rgb(237,224,0, maxColorValue=255),
      rgb(0,30,167, maxColorValue=255),
      rgb(0,208,26, maxColorValue=255),
      rgb(120,0,178, maxColorValue=255),
      rgb(1,235,132, maxColorValue=255),
      rgb(214,0,196, maxColorValue=255),
      rgb(105,180,0, maxColorValue=255),
      rgb(191,88,255, maxColorValue=255),
      rgb(159,255,180, maxColorValue=255),
      rgb(255,3,167, maxColorValue=255),
      rgb(1,244,204, maxColorValue=255),
      rgb(231,41,0, maxColorValue=255),
      rgb(1,220,212, maxColorValue=255),
      rgb(255,90,99, maxColorValue=255),
      rgb(0,129,60, maxColorValue=255),
      rgb(95,93,0, maxColorValue=255),
      rgb(148,0,79, maxColorValue=255),
      rgb(1,124,123, maxColorValue=255),
      rgb(242,188,135, maxColorValue=255),
      rgb(1,200,242, maxColorValue=255),
      rgb(180,80,0, maxColorValue=255),
      rgb(190,203,255, maxColorValue=255),
      rgb(101,15,0, maxColorValue=255))
for(i in 1:length(pca_pop)){
points(tab$EV2[tab$pop == pca_pop[i]],tab$EV1[tab$pop == pca_pop[i]],col = colors[i])
}
legend(0.07,0.03, legend=pca_pop, pch="o", col=colors,pt.cex = 0.7,cex = 0.7)
# dev.off()
# pdf("./global_pair.pdf")
pc.percent = pca$varprop*100
lbls = paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[id_nrow,1:4], col=colors, labels=lbls)
# dev.off()
########################## TREEMIX SECTION ##############################
# Make a plink cluster file matching each individual to a population.
clust = all_fam[all_fam$pop %in% c("San","Yoruba","Mandenka","Papuan","Karitiana","Surui",
                                   "Miaozu","TH-HM","She","CN-JI","Yizu","Naxi",
                                   "Lahu","TH-PL","Bunun","AX-AT","AX-AM","Amis","Paiwan",
                                   "Puyuma","TH-PP","TH-TN","CN-WA","PI-AG","PI-AT",
                                   "PI-AE","PI-IR","PI-MW","PI-MA","PI-UN","PI-UI"),]
clust = clust[,c(1,2,7)]
write.table(clust, file = "./data.clust",row.names = FALSE,col.names = FALSE,quote = FALSE)

system("plink --bfile final_merge_no_LD --freq --missing --within data.clust --out final_merge_no_LD")
system("gzip -c final_merge_no_LD.frq.strat > final_merge_no_LD.frq.strat.gz")
system("python plink2treemix.py final_merge_no_LD.frq.strat.gz final_merge_no_LD.treemix.gz")

## IN SERVER
# treemix -i flip_taitung_combine.treemix.gz -o flip_taitung_combine
# treemix -i flip_taitung_combine.treemix.gz -root San -o flip_taitung_combine

## Small populations TreeMix
