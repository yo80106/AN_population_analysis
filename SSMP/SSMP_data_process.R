setwd("~/SSMP/data")
for(i in 1:22){
    url = paste("http://www.statgen.nus.edu.sg/~SSMP/download/vcf/2012_05/snps/SSM.chr",i,".2012_05.genotypes.vcf.gz",sep="")
    file_name = paste("chr",i,".vcf.gz",sep="")
    print(file_name)
    download.file(url,file_name)
}
system("gunzip *.gz")
for (t in 1:22){
    system(paste("plink -vcf SSM.chr",t,".2012_05.genotypes.vcf --double-id --exclude dot_rm.txt --make-bed --out raw_chr",t,sep = ""))
}
system("plink -vcf SSM.chrX.2012_05.genotypes.vcf --double-id --exclude dot_rm.txt --make-bed --out raw_chr23")
system("plink -vcf SSM.chrY.2013_03.genotypes.vcf --double-id --exclude dot_rm.txt --make-bed --out raw_chr24")
#system("plink -vcf SSM.chrM.2013_04.genotypes.vcf --double-id --make-bed --out chr26")

file_list = list.files("./",pattern = "(.bim|.fam|.bed)$")
file_list = file_list[order(nchar(file_list),file_list)]
t_matrix = matrix(file_list, ncol = 3, byrow = TRUE)
write.table(t_matrix[-1,],"./merge_list.txt",col.names = F,row.names = F,quote = F)

system("plink --bfile raw_chr1 --merge-list merge_list.txt --make-bed --out test")
for(z in 1:nrow(t_matrix)){
    system(paste("plink --bfile raw_chr",z," --exclude test-merge.missnp --make-bed --out chr",z,sep = ""))
}
sub_matrix = gsub("raw_","",t_matrix)
write.table(sub_matrix[-1,], "./sub_merge_list.txt", col.names = F, row.names = F, quote = F)
system("plink --bfile chr1 --merge-list sub_merge_list.txt --make-bed --out ~/SSMP/ssmp")

setwd("~/SSMP")
taitung_bim = read.table("./new_clean_per_snp.bim",stringsAsFactors = FALSE)
ssmp_bim = read.table("./ssmp.bim", stringsAsFactors = FALSE)
t_merge = merge(taitung_bim, ssmp_bim, by = "V2")
write.table(t_merge[,1], "./snp_list.txt", col.names = F, row.names = F, quote = F)

system("plink --bfile new_clean_per_snp --extract snp_list.txt --make-bed --out taitung_extract")
system("plink --bfile ssmp --extract snp_list.txt --make-bed --out ssmp_extract")
system("plink --bfile taitung_extract --bmerge ssmp_extract --make-bed --out test2")
system("plink --bfile taitung_extract --flip test2-merge.missnp --make-bed --out taitung_extract_flip")
system("plink --bfile taitung_extract_flip --bmerge ssmp_extract --make-bed --out test3")
system("plink --bfile taitung_extract_flip --exclude test3-merge.missnp --make-bed --out exclude_taitung_ex_flip")
system("plink --bfile exclude_taitung_ex_flip --bmerge ssmp_extract --recode transpose --out test4")

########################## Malay PCA SECTION ##############################
library(gdsfmt)
library(SNPRelate)

# Conversion from PLINK BED to GDS
bedfile = "./ssmp.bed"
famfile = "./ssmp.fam"
bimfile = "./ssmp.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "ssmp_malay.gds")

genofile = snpgdsOpen("ssmp_malay.gds")
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

pca = snpgdsPCA(genofile)
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
tab = data.frame(sample.id = pca$sample.id,
                 pop = "Malay",
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, col="blue", xlab="PC 2", ylab="PC 1", main="PCA using all SNPs")
legend(0.4,0.1, legend=levels(tab$pop), pch="o", col=1:(nlevels(tab$pop)))
########################## Malay Taiwanese aborigines PCA SECTION ##############################
library(gdsfmt)
library(SNPRelate)

all_fam = read.table("./test4.fam")
all_fam$pop = NA
taitung_list = read.csv("~/PCA/R_analysis/96samples_sex_list.csv",header = F,
                        stringsAsFactors = FALSE)

for(t in 1:nrow(all_fam)){
    for(k in 1:nrow(taitung_list)){
        if(all_fam$V2[t] == taitung_list$V2[k]){
            all_fam$pop[t] = taitung_list$V4[k]
        }
    }
}
all_fam$pop[is.na(all_fam$pop)] = "Malay"


# Conversion from PLINK BED to GDS
bedfile = "./test4.bed"
famfile = "./test4.fam"
bimfile = "./test4.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "ssmp.gds")

genofile = snpgdsOpen("ssmp.gds")
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

pca = snpgdsPCA(genofile)
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
pop = all_fam$pop
tab = data.frame(sample.id = pca$sample.id,
                 pop = factor(pop)[match(pca$sample.id, sample.id)],
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="PC 2", ylab="PC 1", main="PCA using all SNPs")
legend(0.4,0.1, legend=levels(tab$pop), pch="o", col=1:(nlevels(tab$pop)))
########################## SSMP Final merge PCA SECTION ##############################
final_merge_bim = read.table("./final_merge.bim", stringsAsFactors = FALSE)
ssmp_bim = read.table("./ssmp.bim", stringsAsFactors = FALSE)
fm_ssmp_merge = merge(final_merge_bim, ssmp_bim, by = "V2")
write.table(fm_ssmp_merge[,1], "./fm_ssmp_snp_list.txt", col.names = F, row.names = F, quote = F)

system("plink --bfile final_merge --extract fm_ssmp_snp_list.txt --make-bed --out final_merge_extract")
system("plink --bfile ssmp --extract fm_ssmp_snp_list.txt --make-bed --out fm_ssmp_extract")
system("plink --bfile final_merge_extract --bmerge fm_ssmp_extract --make-bed --out fm_ssmp_merge_1")
system("plink --bfile final_merge_extract --flip fm_ssmp_merge_1-merge.missnp --make-bed --out fm_ssmp_merge_2")
system("plink --bfile fm_ssmp_merge_2 --bmerge fm_ssmp_extract --make-bed --out fm_ssmp_merge_3")

library(gdsfmt)
library(SNPRelate)
library(stringr)

all_fam = read.table("./fm_ssmp_merge_3.fam")
all_fam$pop = NA
taitung_list = read.csv("~/PCA/R_analysis/96samples_sex_list.csv",header = F,
                        stringsAsFactors = FALSE)

# Fill in Pan-Asia population list
for(i in 1:nrow(all_fam)){
    pop =  str_sub(all_fam$V2[i],start = 1,end = 5)
    all_fam$pop[i] = pop 
}

# Fill in Taitung population list
for(t in 1:nrow(all_fam)){
    for(k in 1:nrow(taitung_list)){
        if(all_fam$V2[t] == taitung_list$V2[k]){
            all_fam$pop[t] = taitung_list$V4[k]
        }
    }
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
unique(all_fam$pop)
# Fill in SSMP population list
all_fam$pop[grep("SSM", all_fam$pop)] = "Malay"


# Conversion from PLINK BED to GDS
bedfile = "./fm_ssmp_merge_3.bed"
famfile = "./fm_ssmp_merge_3.fam"
bimfile = "./fm_ssmp_merge_3.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "fm_ssmp.gds")

genofile = snpgdsOpen("fm_ssmp.gds")
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

pca = snpgdsPCA(genofile)
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
pop = all_fam$pop
tab = data.frame(sample.id = pca$sample.id,
                 pop = factor(pop)[match(pca$sample.id, sample.id)],
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)
pdf("./ssmp_final_merge.pdf")
plot(tab$EV2, tab$EV1,col = rgb(214,214,214, maxColorValue = 255), xlab="PC 2", ylab="PC 1", main="PCA using all SNPs")
pca_pop = c("Paiwan","Amis","Bunun","Puyuma","TW-HA",
            "TW-HB","SG-CH","SG-ML", "SG-ID","Malay",
            "MY-KS","MY-JH","MY-TM","MY-KN","MY-MN",
            "PI-AG","PI-AT","Papuan","NANMelanesian")
id_nrow = which(tab$pop %in% pca_pop)
id_row_pop = tab$pop[tab$pop %in% pca_pop]
colors = c("gold2","darkorange","darkgoldenrod3",
          "yellowgreen","green4","antiquewhite1",
          "lavenderblush3","sienna3","tomato4",
          "palevioletred1","cyan3","red2",
          "navy","pink1","purple","slategray1",
          "cadetblue","royalblue1","violetred1")
for(i in 1:length(pca_pop)){
    points(tab$EV2[tab$pop == pca_pop[i]],tab$EV1[tab$pop == pca_pop[i]], col = colors[i],lwd = 3)
}
legend(-0.085,0, legend=c(pca_pop[-19],"Melanesian"), pch=1, col=colors,pt.cex = 0.8,cex = 0.8, ncol = 3, pt.lwd = 3)
dev.off()

pdf("./ssmp_final_merge2.pdf")
plot(tab$EV2, tab$EV1,col = rgb(214,214,214, maxColorValue = 255), xlab="PC 2", ylab="PC 1", main="PCA using all SNPs")
pca_pop = c("TW-HA","TW-HB","SG-CH","SG-ML","SG-ID",
            "Malay","MY-KS","MY-JH","PI-AG","PI-AT",
            "Papuan","NANMelanesian")
id_nrow = which(tab$pop %in% pca_pop)
id_row_pop = tab$pop[tab$pop %in% pca_pop]
colors = c("gold2","darkorange","darkgoldenrod3",
           "yellowgreen","green4","antiquewhite1",
           "cadetblue","sienna3","tomato4",
           "palevioletred1","cyan3","red2")
for(i in 1:length(pca_pop)){
    points(tab$EV2[tab$pop == pca_pop[i]],tab$EV1[tab$pop == pca_pop[i]], col = colors[i],lwd = 3)
}
legend(-0.085,0, legend=c(pca_pop[-12],"Melanesian"), pch=1, col=colors,pt.cex = 0.8,cex = 0.8, ncol = 3, pt.lwd = 3)
dev.off()

# Make a plink cluster file matching each individual to a population.
clust = all_fam[,c(1,2,7)]
write.table(clust, file = "./data.clust",row.names = FALSE,col.names = FALSE,quote = FALSE)

system("plink --bfile fm_ssmp_merge_3 --freq --missing --within data.clust --out fm_ssmp_merge_3")
system("gzip -c fm_ssmp_merge_3.frq.strat > fm_ssmp_merge_3.frq.strat.gz")
system("python plink2treemix.py fm_ssmp_merge_3.frq.strat.gz fm_ssmp_merge_3.treemix.gz")

## IN SERVER
# treemix -i fm_ssmp_merge_3.treemix.gz -root San -o fm_ssmp_merge_3 

## Venn diagram
library(VennDiagram)
nssmp_bim = nrow(ssmp_bim)
nfinal_bim = nrow(final_merge_bim)
a = merge(ssmp_bim,final_merge_bim,by = "V2")
n12 = length(unique(a$V2))

pdf("./venn.pdf")
draw.pairwise.venn(area1 = nssmp_bim, area2 = nfinal_bim, n12, fill = c("skyblue", "pink1"), cex = 2)
dev.off()

########################## SSMP SG MY PCA SECTION ##############################
malay_fam = all_fam[(all_fam$pop %in% c("Malay","SG-CH","SG-ML","SG-ID",
                                        "Paiwan","Amis","Bunun","Puyuma",
                                        "AX-AM","AX-AT","TW-HA","TW-HB",
                                        "MY-TM","MY-BD","MY-KS","MY-JH",
                                        "MY-MN","MY-KN","PI-AG","PI-AT")),]
write.table(malay_fam[,c(1,2)],file = "./malay_sg_pop_list.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)
system("plink --bfile fm_ssmp_merge_3 --keep malay_sg_pop_list.txt --make-bed --out malay_sg")

library(gdsfmt)
library(SNPRelate)

fam = read.table("./malay_sg.fam")
ind_info = malay_fam
sum(as.character(ind_info$V2)!=as.character(fam$V2))

bedfile = "./malay_sg.bed"
famfile = "./malay_sg.fam"
bimfile = "./malay_sg.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "malay_sg.gds")

genofile = snpgdsOpen("malay_sg.gds")
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

pca = snpgdsPCA(genofile)
sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
pop = malay_fam$pop
tab = data.frame(sample.id = pca$sample.id,
                 pop = factor(pop)[match(pca$sample.id, sample.id)],
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)
palette(c(V1="gold2",V2="darkorange",V3="darkgoldenrod3",
          V4="yellowgreen",V5="green4",v6="antiquewhite1",
          V7="lavenderblush3",V8="sienna3",V9="tomato4",
          V10="palevioletred1",V11="cyan3",V12="red2",
          V13="navy",V14="pink1",V15="purple",
          V16="slategray1",V17="cadetblue",V18="royalblue1",
          V19="violetred1",V20="coral"))
        #V21="chartreuse1",V22="darkorchid4",V23="olivedrab4",V24="aquamarine1"
pdf("./pca_malay_sg2.pdf")
plot(tab$EV1, tab$EV2, col=tab$pop, xlab="PC 1", ylab="PC 2", main="PCA using all SNPs", lwd = 3)
# text(tab$EV1,tab$EV2,labels = tab$sample.id,cex = 0.4)
legend(-0.05,0.195, legend=levels(tab$pop), pch=1, col=1:(nlevels(tab$pop)),pt.cex = 0.8,cex = 0.8, ncol = 4, pt.lwd = 3)
dev.off()
#for(n in 1:length(file_list)){
    temp_bim = read.table(file_list[n],header = F,stringsAsFactors = FALSE)
    temp_taitung = subset(taitung_bim, V1 == n)
    for (n2 in 1:nrow(temp_bim)) {
        for (n3 in 1:nrow(taitung_bim)) {
            if(temp_bim$V4[n2] %in% temp_taitung$V4[n3]){
                temp_bim$V2[n2] = temp_taitung$V2[n3]
            }
        }
    }
}
#for (n in 1:length(file_list)) {
    if(!exists("ssmp_map")){
        ssmp_map = read.table(file_list[n],header = F,stringsAsFactors = FALSE)
        ssmp_map = ssmp_map[ssmp_map$V2 != ".",]
    }
    if(exists("ssmp_map")){
        temp_ssmp = read.table(file_list[n],header = F,stringsAsFactors = FALSE)
        temp_ssmp = temp_ssmp[temp_ssmp$V2 != ".",]
        ssmp_map = rbind(ssmp_map,temp_ssmp)
        rm(temp_ssmp)
    }
}


