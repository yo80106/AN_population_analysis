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

missing_het_ind = function(input.name, pop.list, output.name){
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
    
    sex = read.table(paste0(pop.list,".txt"), header = F)
    pop = names(table(sex[,3]))
    pop_col = rainbow(length(pop))
    colors = densCols(imiss$logF_MISS,het$Het_propo)
    for(i in 1:length(sex[,2])){
        if(sex[,2][i] == het$IID[i]){
            colors[i] = pop_col[sex[,3][i]]
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
    system(paste0("plink --bfile ",input.name," --indep-pairwise 50 5 0.2 --out ",output.name))
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

ind_qc_info = function(input.name, output.name){
    # if(sum(installed.packages(.libPaths())[,1] %in% "venneuler") == 0){
    #     install.packages("venneuler")
    # }
    # library(venneuler)
    problem_all = c()
    if(length(list.files(pattern = paste0(input.name,"_sex_problem.list")) == 1)){
        ind_sex = read.table(paste0(input.name, "_sex_problem.list"), stringsAsFactors = F)
        print("Loading sex problem list.")
        problem_all = c(problem_all, Sex = ind_sex)
    }else if(length(list.files(pattern = paste0(input.name,"_sex_problem.list")) == 0)){
        print("No sex problem list!")
    }else{
        print("There are something wrong with sex problem list.")
    }
    if(length(list.files(pattern = paste0(input.name,"_ibd_problem.list")) == 1)){
        ind_ibd = read.table(paste0(input.name, "_ibd_problem.list"), stringsAsFactors = F)
        print("Loading ibd problem list.")
        problem_all = c(problem_all, IBD = ind_ibd)
    }else if(length(list.files(pattern = paste0(input.name,"_sex_problem.list")) == 0)){
        print("No ibd problem list!")
    }else{
        print("There are something wrong with ibd problem list.")
    }
    if(length(list.files(pattern = paste0(input.name,"_miss_het_problem.list")) == 1)){
        ind_miss = read.table(paste0(input.name, "_miss_het_problem.list"), stringsAsFactors = F)
        print("Loading het and miss problem list.")
        problem_all = c(problem_all, Miss = ind_miss)
    }else if(length(list.files(pattern = paste0(input.name,"_miss_het_problem.list")) == 0)){
        print("No het and miss problem list!")
    }else{
        print("There are something wrong with miss problem list.")
    }
    if(length(problem_all) == 6){
        ABC = Reduce(intersect, c(problem_all[1], problem_all[3], problem_all[5]))
        AB = Reduce(intersect, c(problem_all[1], problem_all[3]))
        AC = Reduce(intersect, c(problem_all[1], problem_all[5]))
        BC = Reduce(intersect, c(problem_all[3], problem_all[5]))
        out_table = list(ABC, AB, AC, BC)
        
        A = sub("\\..*", "", names(problem_all[1]))
        B = sub("\\..*", "", names(problem_all[3]))
        C = sub("\\..*", "", names(problem_all[5]))
        names(out_table) = c(paste(A, B, C, sep = "+"), paste(A, B, sep = "+"), paste(A, C, sep = "+"), paste(B, C, sep = "+"))
        df = List_to_DF(out_table)
        write.csv(df, file = paste0(output.name, "_ind_qc_info.csv"), quote = F, row.names = F)
        # AB_num = length(out_table[[2]])
        # AC_num = length(out_table[[3]])
        # BC_num = length(out_table[[4]])
        # ABC_num = length(out_table[[1]])
        # A_num = length(problem_all[[1]]) - AB_num - AC_num + ABC_num
        # B_num = length(problem_all[[3]]) - AB_num - BC_num + ABC_num
        # C_num = length(problem_all[[5]]) - AC_num - BC_num + ABC_num
        # threeVen = venneuler(c(A = A_num , B = B_num, C = C_num,
        #                   "A&B" = AB_num, "A&C" = AC_num, 
        #                      "B&C" = BC_num, "A&B&C" = ABC_num))
        # threeVen$labels = c(
        #     paste0(A, "\n", length(problem_all[[1]])),
        #     paste0(B, "\n", length(problem_all[[3]])),
        #     paste0(C, "\n", length(problem_all[[5]]))
        # )
        # pdf(paste0(output.name, "_threeVen.pdf"))
        # plot(threeVen)
        # dev.off()
    }else if(length(problem_all) == 4){
        AB = Reduce(intersect, c(problem_all[1], problem_all[3]))
        A = sub("\\..*", "", names(problem_all[1]))
        B = sub("\\..*", "", names(problem_all[3]))
        df = data.frame(AB)
        colnames(df) = paste(A, B, sep = "+")
        write.csv(df, file = paste0(output.name, "_ind_qc_info.csv"), quote = F, row.names = F)
        # twoVen = venneuler(c(A = length(problem_all[[1]]), B = length(problem_all[[3]]), "A&B" = length(out_table[[2]])))
        # twoVen$labels = c(
        #     paste0(A, "\n", length(problem_all[[1]])),
        #     paste0(B, "\n", length(problem_all[[3]]))
        # )
        # pdf(paste0(output.name, "_tewVen.pdf"))
        # plot(twoVen)
        # dev.off()
    }
}

List_to_DF = function(input.list){
    lengths = lapply(input.list, length)
    max_num = as.integer(lengths[which.max(lengths)])
    new_df = matrix(nrow = max_num, ncol = length(input.list))
    for(n in 1:length(input.list)){
        new_df[,n] = c(unlist(input.list[n]), rep("-", max_num-lengths(input.list[n])))
    }
    new_df = as.data.frame(new_df, stringsAsFactors = F)
    colnames(new_df) = c(names(input.list))
    return(new_df)
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