my_pca = function(plink.file, data.clust, pop.list, output.name){
    if(sum(installed.packages(.libPaths())[,1] %in% "SNPRelate") == 0){
        source("https://bioconductor.org/biocLite.R")
        biocLite("SNPRelate")
    }else if(sum(installed.packages(.libPaths())[,1] %in% "randomcoloR") == 0){
        install.packages("randomcoloR")
    }else if(sum(installed.packages(.libPaths())[,1] %in% "gdsfmt") == 0){
        install.packages("gdsfmt")
    }
    library(gdsfmt)
    library(SNPRelate)
    library(randomcoloR)
    all_ind = read.table(data.clust, stringsAsFactors = FALSE)
    select_pop = read.table(pop.list, stringsAsFactors = FALSE)
    select_ind = all_ind[all_ind$V3 %in% select_pop$V1,]
    write.table(select_ind[,c(1,2)], file = paste0(output.name, ".ind.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    system(paste0("plink --bfile ", plink.file, " --keep ",output.name,".ind.txt", " --make-bed --out ", output.name))
    
    bedfile = paste0(output.name, ".bed")
    famfile = paste0(output.name, ".fam")
    bimfile = paste0(output.name, ".bim")
    snpgdsBED2GDS(bedfile, famfile, bimfile, paste0(output.name, ".gds"))
    
    genofile = snpgdsOpen(paste0(output.name, ".gds"))
    pca = snpgdsPCA(genofile)
    sample.id = read.gdsn(index.gdsn(genofile, "sample.id"))
    pop = select_ind$V3
    pc.percent = round(pca$varprop*100, 2)
    tab = data.frame(sample.id = pca$sample.id,
                     pop = factor(pop)[match(pca$sample.id, sample.id)],
                     EV1 = pca$eigenvect[,1],    # the first eigenvector
                     EV2 = pca$eigenvect[,2],   # the second eigenvector
                     varprop = pc.percent
    )
    write.table(tab, file = paste0(output.name, "_pca.tab"), quote = F, row.names = F, col.names = F)
    palette(distinctColorPalette(k = nlevels(tab$pop)))
    pdf(file = paste0(output.name, "_pca.pdf"))
    par(xpd = T, mar = par()$mar + c(0,0,0,5))
    plot(tab$EV2, tab$EV1, col=tab$pop, xlab="PC 2", ylab="PC 1", lwd = 3)
    legend("right", legend = levels(tab$pop), col = 1:(nlevels(tab$pop)), inset=c(-0.2,0.2), pch = 1, pt.cex = 0.8, cex = 0.8, pt.lwd = 3)
    dev.off()
    pdf(file = paste0(output.name, "_pairs.pdf"))
    pc.percent = pca$varprop*100
    lbls = paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
    pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)
    dev.off()
    pdf(file = paste0(output.name, "_var.pdf"))
    plot(tab$varprop[complete.cases(tab$varprop)], type = "b" , ylab = "Variance(%)", xlab = "PCs")
    dev.off()
    
}
