########################## DATA PREPROCESSING ##############################
library(stringr)
hgdp_pan_fam = read.table("~/PCA/R_analysis/HGDP-Pan-LD.fam")

hgdp_pan_fam["pop"] = NA
for(i in 1:nrow(hgdp_pan_fam)){
   pop =  str_sub(hgdp_pan_fam$V1[i],start = 1,end = 5)
    hgdp_pan_fam$pop[i] = pop 
}

hgdp_list = read.table("~/PCA/R_analysis/HGDP-country-list.txt",header = T,stringsAsFactors = F)
parse = unlist(strsplit(hgdp_list$ID_Pop,"_"))
hgdp_list["HGDP_ID"] = parse[seq(1,length(parse),2)]
hgdp_list["pop"] = parse[seq(2,length(parse),2)]

for(t in 1:nrow(hgdp_pan_fam)){
    for(k in 1:nrow(hgdp_list)){
        if(hgdp_pan_fam$V1[t] == hgdp_list$HGDP_ID[k]){
            hgdp_pan_fam$pop[t] = hgdp_list$pop[k]
        }
    }
}

########################## PCA SECTION ##############################
source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")

library(gdsfmt)
library(SNPRelate)

philippines = c("AX-AM","AX-AT","Lahu","Naxi","TH-HM","TH-PP",
                "TH-TN","Papuan","PI-AT",
                "PI-MW","PI-MA","PI-UN","PI-UI","PI-UB")

ind_info_num = unlist(lapply(philippines, grep, hgdp_pan_fam$pop))
ind_info = hgdp_pan_fam[ind_info_num,]
table(ind_info$pop)

write.table(ind_info[,c(1,2)],file = "./philippines_ind.txt",quote = FALSE,
            row.names = FALSE,col.names = FALSE)

# Conversion from PLINK BED to GDS
bedfile = "~/PCA/R_analysis/philippines.2.bed"
famfile = "~/PCA/R_analysis/philippines.2.fam"
bimfile = "~/PCA/R_analysis/philippines.2.bim"
snpgdsBED2GDS(bedfile, famfile, bimfile, "IPCgeno.gds")

genofile = snpgdsOpen("IPCgeno.gds")
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
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), xlab="PC 1", ylab="PC 2", 
     pch = as.integer(tab$pop),main="PCA using all SNPs")
legend("bottomright", legend=levels(tab$pop), 
       pch=1:(nlevels(tab$pop)), col=1:(nlevels(tab$pop)))

pc.percent = pca$varprop*100
lbls = paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)
