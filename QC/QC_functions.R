sex_check = function(input.name, output.name){
    system(paste0("plink --bfile ",input.name," --check-sex --out ",output.name))
    sexcheck = read.table(paste0(output.name,".sexcheck"), header = T)
    colors = densCols(sexcheck$F)
    pdf("sex_distribution.pdf")
    par(cex.lab=1.5,cex.main=2.5,cex.axis =1.7)
    plot(sexcheck$PEDSEX,sexcheck$F,pch=19, xlim = c(0,2),
         col=colors,axes = F, xlab = "Sex record", ylab = "Homozygote rate")
    axis(1,at=seq(0,2,1),labels = c("Undetermined","Male","Female"))
    axis(2,at=seq(0,1,0.1))
    abline(h=0.8,lty=2,col="green",lwd=3)
    abline(h=0.2,lty=2,col="green",lwd=3)
    title("Distribution of individual's sex")
    dev.off()
    fail_sex = which(sexcheck$STATUS == "PROBLEM")
    fail_sex_id = sexcheck[c(fail_sex),c(1,2)]
    write.table(fail_sex_id, paste0(output.name, "_sex_problem.list"), col.names = F, row.names = F, quote = F)
}

missing_het_ind = function(input.name, sex.list, output.name){
    system(paste0("plink --bfile ",input.name," --missing --out ",output.name))
    system(paste0("plink --bfile ",input.name," --het --out ",output.name))
    
    imiss = read.table(paste0(output.name, ".imiss"), header = T)
    het = read.table(paste0(output.name, ".het"), header = T)
    
    # Calculate call rate, and add a column
    imiss$CALL_RATE = 1-imiss$F_MISS 
    # Log 10 the F_MISS column, add a column
    imiss$logF_MISS = log10(imiss[,6]) 
    # Calculate heterozygosity rate, and add a column
    het$Het_propo = (het$N.NM. - het$O.HOM.)/het$N.NM.
    # Find out 'NaN' value in meanHet column
    # het$Het_propo = ifelse(het$Het_propo=="NaN", 0,het$Het_propo) 
    
    read.csv(sex.list, header = F)
    pop = names(table(sex$V4))
    pop_col = rainbow(length(pop))
    colors = densCols(imiss$logF_MISS,het$Het_propo)
    for(i in 1:length(sex$V2)){
        if(sex$V2[i] == het$IID[i]){
            colors[i] = pop_col[sex$V4[i]]
        }
    }
    pdf("imiss-vs-het.pdf")
    plot(imiss$logF_MISS,het$Het_propo, col=colors,pch=20, xlim=c(-3,0),
         ylim=c(0,1),xlab="Proportion of missing genotypes", 
         ylab="Heterozygosity rate",axes=FALSE)
    axis(2,at=seq(0,1,0.1),tick=T)
    axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
    legend("bottomright", legend=pop, col=pop_col, pch=16, ncol=2, bty="n", cex=0.8)
    # Heterozygosity thresholds (Horizontal Line) +-3 s.d from the mean
    abline(h=mean(het$Het_propo)-(3*sd(het$Het_propo)),col="RED",lty=2)
    abline(h=mean(het$Het_propo)+(3*sd(het$Het_propo)),col="RED",lty=2)
    # Missing Data Thresholds (Vertical Line)
    abline(v=log10(0.03), col="BLUE", lty=2) # THRESHOLD=0.03
    title("Distribution of missingness\nand heterozygosity scores")
    dev.off()
    # Use R to identify individuals with high missingness and/or outlier 
    # heterozygosity based on preselected cutoff
    imiss_het = merge(het,imiss,by = "FID")
    fail_imisshet = imiss_het$Het_propo < mean(het$Het_propo)-(3*sd(het$Het_propo)) | imiss_het$Het_propo > mean(het$Het_propo)+(3*sd(het$Het_propo)) | imiss_het$F_MISS >= 0.03
    which_bad = which(fail_imisshet,TRUE)
    fail_imisshet_id = imiss_het[c(which_bad),c(1,2)]
    colnames(fail_imisshet_id) = c("FID","IID")
    write.table(fail_imisshet_id,file = paste0(output.name, "_miss_het_problem.list"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}

IBD = function(input.name, output.name){
    system(paste0("plink --file ",input.name," --indep-pairwise 50 5 0.2 --out ",output.name))
    system(paste0("plink --bfile ",input.name," --extract ",output.name,".prune.in --genome --out ",output.name))
    system(paste0("plink --bfile ",input.name," --missing --out ",output.name))
    ibd_data = read.table(paste0(output.name, ".genome"), header = T)
    
    # If there are any 'NaN' in data.
    # bad = which(is.nan(ibd_data$PI_HAT))
    # ibd_new = ibd_data$PI_HAT[-bad]
    pro = data.frame()
    total = length(ibd_data$PI_HAT)
    for(i in 1:length(unique(ibd_data$PI_HAT))){
        pro[i,1] = unique(ibd_data$PI_HAT)[i]
        pro[i,2] = sum(ibd_data$PI_HAT == unique(ibd_data$PI_HAT)[i])
        pro[i,3] = sum(ibd_data$PI_HAT == unique(ibd_data$PI_HAT)[i])/total
    }
    colnames(pro) = c("ibd","total","proportion")
    pdf("IBD.pdf")
    plot(pro$ibd, pro$proportion, type="h", ylim=c(0,0.8), 
         xlim=c(0,1),ylab="probability",xlab="IBD", lwd=8, col="orange")
    abline(v=0.1875, lty=2, col="red")
    text(x=0.2280, y=0.4, labels=paste("cutoff=",0.1875, sep=""), col="red")
    title("Propotion of the different IBD")
    dev.off()
    # Use R to identify individuals with high IBD
    imiss_data = read.table(paste0(output.name, ".imiss"), header = T)
    first = ibd_data[,1][which(ibd_data$PI_HAT > 0.1875)]
    second = ibd_data[,3][which(ibd_data$PI_HAT > 0.1875)]
    
    fail_ibd = c()
    for(n in 1:length(first)){
        if(imiss_data[first[n],6] > imiss_data[second[n],6]){
            fail_ibd[n] = imiss_data[first[n],1]
        }else if(imiss_data[first[n],6] < imiss_data[second[n],6]){
            fail_ibd[n] = imiss_data[second[n],1]
        }else if(imiss_data[first[n],6] == imiss_data[second[n],6]){
            fail_ibd[n] = imiss_data[first[n],1]
        }
    }
    fail_ibd = sort(unique(fail_ibd))
    fail_ibd_id = imiss_data[c(fail_ibd),c(1,2)]
    write.table(fail_ibd_id, file = paste0(output.name, "_ibd_problem.list"), col.names = F, row.names = F, quote = F)
}

ind_qc_rm = function(input.name, output.name){
    problem_lists = list.files("./", pattern = paste0("^",input.name,"(.*)problem.list$"))
    if(length(problem_lists == 3)){
        a = c()
        for(i in 1:3){
            b = read.table(problem_lists[i], stringsAsFactors = F)
            a = rbind(a,b)
        }
    }else{
        print("There are more than three problem list files in your directory!")
        }
    fail_qc_inds = write.table(unique(a[order(a$V1),]), paste0(output.name, "_fail_ind_qc.txt"), row.names = F, col.names = F, quote = F)   
}

missing_snp = function(input.name, output.name){
    system(paste0("plink --bfile ",input.name," --missing --out ",output.name))
    # Plot snp missing distribution
    # Load SNP frequency file and generate histogram
    b.frq = read.table(paste0(output.name, ".lmiss"), header = T)
    pdf("snpmiss_plot.pdf")
    plot(ecdf(b.frq$F_MISS), xlim = c(0,0.10), ylim = c(0,1), pch = 20, 
         main = "SNP Missingness Distribution", xlab = "Missingness Frequency", 
         ylab = "Fraction of SNPs", col = "blue")
    abline(v=0.03, col = "red", lty = 2)
    dev.off()
}

hwe_test = function(input.name, output.name){
    system(paste0("plink --bfile ",input.name," --hardy --out ",output.name))
    hwe.file = read.table(paste0(output.name, ".hwe"), header = T)
    pdf("hwe_p_value.pdf")
    plot(ecdf(na.omit(hwe.file$P)), xlim = range(hwe.file$P,na.rm = T), 
         ylim = c(0,1), pch = 20, main = "Distribution of HWE P-Value", 
         xlab = "P-value", ylab = "Fraction of SNPs",col="blue")
    abline(v = 0.00001, col = "red", lty = 2)
    dev.off()
}