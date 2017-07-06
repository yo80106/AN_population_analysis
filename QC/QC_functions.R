sex_check = function(input.name, female, male){
    if(!missing(male) & !missing(female)){
        system(paste0("plink --bfile ",input.name," --check-sex ",female, " ", male, " --out ", input.name))
    }else if(missing(male) & missing(female)){
        system(paste0("plink --bfile ",input.name," --check-sex --out ", input.name))
    }else{
        print("You should either give the threholds for female and male repectively or leave both blanks.")
    }
sexcheck = read.table(paste0(input.name,".sexcheck"), header = T)
colors = densCols(sexcheck$F)
pdf("sex_distribution.pdf")
par(cex.lab=1.5,cex.main=2.5,cex.axis =1.7)
plot(sexcheck$PEDSEX,sexcheck$F,pch=19, xlim = c(0,2),
     col=colors,axes = F, xlab = "Sex record", ylab = "Homozygote rate")
axis(1,at=seq(0,2,1),labels = c("Undetermined","Male","Female"))
axis(2,at=seq(0,1,0.1))
abline(h=0.8,lty=2,col="green",lwd=3)
abline(h=0.25,lty=2,col="green",lwd=3)
title("Distribution of individual's sex")
dev.off()
fail_sex = which(sexcheck$STATUS == "PROBLEM")
fail_sex_id = sexcheck[c(fail_sex),c(1,2)]
write.table(fail_sex_id, paste0(input.name, "_sex_problem.list"), col.names = F, row.names = F, quote = F)
}

missing_het_ind = function(input.name, pop.list, missing_rate){
    system(paste0("plink --bfile ",input.name," --missing --out ", input.name))
    system(paste0("plink --bfile ",input.name," --het --out ", input.name))
    
    imiss = read.table(paste0(input.name, ".imiss"), header = T)
    het = read.table(paste0(input.name, ".het"), header = T)
    
    # Calculate call rate, and add a column
    imiss$CALL_RATE = 1-imiss$F_MISS 
    # Log 10 the F_MISS column, add a column
    imiss$logF_MISS = log10(imiss[,6]) 
    # Calculate heterozygosity rate, and add a column
    het$Het_propo = (het$N.NM. - het$O.HOM.)/het$N.NM.
    # Find out 'NaN' value in meanHet column
    # het$Het_propo = ifelse(het$Het_propo=="NaN", 0,het$Het_propo) 
    if(!missing(pop.list)){
        sex = read.table(paste0(pop.list), header = F)
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
        abline(v=log10(missing_rate), col="BLUE", lty=2) # THRESHOLD=missing_rate
        title("Distribution of missingness\nand heterozygosity scores")
        dev.off()
    }else if(missing(pop.list)){
        colors = densCols(imiss$logF_MISS,het$Het_propo)
        pdf("imiss-vs-het.pdf")
        plot(imiss$logF_MISS,het$Het_propo, col=colors,pch=20, xlim=c(-3,0),
             ylim=c(0,1),xlab="Proportion of missing genotypes", 
             ylab="Heterozygosity rate",axes=FALSE)
        axis(2,at=seq(0,1,0.1),tick=T)
        axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
        # Heterozygosity thresholds (Horizontal Line) +-3 s.d from the mean
        abline(h=mean(het$Het_propo)-(3*sd(het$Het_propo)),col="RED",lty=2)
        abline(h=mean(het$Het_propo)+(3*sd(het$Het_propo)),col="RED",lty=2)
        # Missing Data Thresholds (Vertical Line)
        abline(v=log10(missing_rate), col="BLUE", lty=2) # THRESHOLD=missing_rate
        title("Distribution of missingness\nand heterozygosity scores")
        dev.off()
    }
    # Use R to identify individuals with high missingness and/or outlier 
    # heterozygosity based on preselected cutoff
    imiss_het = merge(het,imiss,by = "FID")
    fail_imisshet = imiss_het$Het_propo < mean(het$Het_propo)-(3*sd(het$Het_propo)) | imiss_het$Het_propo > mean(het$Het_propo)+(3*sd(het$Het_propo)) | imiss_het$F_MISS >= missing_rate
    which_bad = which(fail_imisshet,TRUE)
    fail_imisshet_id = imiss_het[c(which_bad),c(1,2)]
    colnames(fail_imisshet_id) = c("FID","IID")
    write.table(fail_imisshet_id,file = paste0(input.name, "_miss_het_problem.list"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}

IBD = function(input.name, ibd, mut_error){
    system(paste0("plink --bfile ", input.name," --missing --out ", input.name))
    ibd_data = read.table(paste0(input.name, ".genome"), header = T)
    imiss_data = read.table(paste0(input.name, ".imiss"), header = T)
    
    # If there are any 'NaN' in data.
    # bad = which(is.nan(ibd_data$PI_HAT))
    # ibd_new = ibd_data$PI_HAT[-bad]
    par(mfrow = c(1,1))
    pdf(paste0(input.name, "_IBD.pdf"))
    hist_info = hist(ibd_data$PI_HAT, plot = F)
    hist_info$prob = hist_info$counts/sum(hist_info$counts)
    plot(hist_info$mids, hist_info$prob, xlim = c(0,1), ylim = c(0,1), type = "h", col = "blue", 
         xlab = "IBD", ylab = "Proportion of IBD", main = paste0(input.name, " n=", nrow(imiss_data)))
    dev.off()
    
    IBD_inner = function(t.ibd, imiss, v.ibd){
        first = t.ibd[,1][which(t.ibd$PI_HAT >= v.ibd)]
        second = t.ibd[,3][which(t.ibd$PI_HAT >= v.ibd)]
        fail_ibd = c()
        if(length(first) == 0){
            message = paste0("There are no failed individuals under ", v.ibd, " cutoff.")
            cat(message, file= paste0(input.name, "_no_fail_ibd.txt"))
        }else{
            for(i in 1:length(first)){
                if(imiss[which(imiss[,1] == first[i]),6] > imiss[which(imiss[,1] == second[i]),6]){
                    fail_ibd[i] = first[i]
                }else if(imiss[which(imiss[,1] == first[i]),6] < imiss[which(imiss[,1] == second[i]),6]){
                    fail_ibd[i] = second[i]
                }else if(imiss[which(imiss[,1] == first[i]),6] == imiss[which(imiss[,1] == second[i]),6]){
                    fail_ibd[i] = first[i]
                }
            }
            fail_ibd = sort(unique(fail_ibd))
            fail_ibd_id = imiss[which(imiss[,1] %in% c(fail_ibd)),c(1,2)]
            write.table(fail_ibd_id, file = paste0(input.name, "_fail_ibd.list"), col.names = F, row.names = F, quote = F)   
        }
    }
    # pro = data.frame()
    # total = length(ibd_data$PI_HAT)
    # for(i in 1:length(unique(ibd_data$PI_HAT))){
    #     pro[i,1] = unique(ibd_data$PI_HAT)[i]
    #     pro[i,2] = sum(ibd_data$PI_HAT == unique(ibd_data$PI_HAT)[i])
    #     pro[i,3] = sum(ibd_data$PI_HAT == unique(ibd_data$PI_HAT)[i])/total
    # }
    # colnames(pro) = c("ibd","total","proportion")
    # pdf("IBD.pdf")
    # plot(pro$ibd, pro$proportion, type="h", ylim=c(0,0.8), 
    #      xlim=c(0,1),ylab="probability",xlab="IBD", lwd=8, col="orange")
    # abline(v=ibd, lty=2, col="red")
    # text(x=0.2280, y=0.4, labels=paste("cutoff=",ibd, sep=""), col="red")
    # title("Propotion of the different IBD")
    # dev.off()
    # Use R to identify individuals with high IBD
    if(missing(ibd)){
        if(missing(mut_error)){
            if(quantile(ibd_data$PI_HAT, 0.95)[[1]] < 0.25){
                IBD_inner(t.ibd = ibd_data, imiss = imiss_data, v.ibd = 0.25)
                txt = paste0("95% quantile value of IBD in ", input.name, " is < 0.25, so the cutoff was set to 0.25.")
                cat(txt, file= paste0(input.name, "_IBD_log.txt"))
            }else if(quantile(ibd_data$PI_HAT, 0.95)[[1]] >= 0.25 & quantile(ibd_data$PI_HAT, 0.95)[[1]] < 0.5){
                IBD_inner(t.ibd = ibd_data, imiss = imiss_data, v.ibd = 0.25)
                txt = paste0("95% quantile value of IBD in ", input.name, " is < 0.5, so the cutoff was set to 0.5.")
                cat(txt, file= paste0(input.name, "_IBD_log.txt"))
            }
        }else if(!missing(mut_error)){
            if(quantile(ibd_data$PI_HAT, 0.95)[[1]] < 0.25-mut_error){
                IBD_inner(t.ibd = ibd_data, imiss = imiss_data, v.ibd =  0.25-mut_error)
                txt = paste0("95% quantile value of IBD in ", input.name, " is < 0.25, so the cutoff was set to ",  0.25-mut_error, ".")
                cat(txt, file= paste0(input.name, "_IBD_log.txt"))
            }else if(quantile(ibd_data$PI_HAT, 0.95)[[1]] >=  0.25-mut_error & quantile(ibd_data$PI_HAT, 0.95)[[1]] < 0.5-mut_error){
                IBD_inner(t.ibd = ibd_data, imiss = imiss_data, v.ibd = 0.5-mut_error)
                txt = paste0("95% quantile value of IBD in ", input.name, " is < 0.5, so the cutoff was set to ", 0.5-mut_error, ".")
                cat(txt, file= paste0(input.name, "_IBD_log.txt"))
            }
        }
    }else if(!missing(ibd)){
        IBD_inner(t.ibd = ibd_data, imiss = imiss_data, v.ibd = ibd)
    }
}

IBD_rm = function(output.name){
    problem_lists = list.files("./", pattern = "_fail_ibd.list$")
    a = c()
    for(i in 1:length(problem_lists)){
        b = read.table(problem_lists[i], stringsAsFactors = F)
        print(problem_lists[i])
        a = rbind(a,b)
    }
    write.table(a, file = paste0(output.name, "_ibd_problem.list"), row.names = F, col.names = F, quote = F)
}

IBD_pop = function(input.name, cutoff, data.clust){
    system(paste0("plink --bfile ",input.name," --indep-pairwise 50 5 0.2 --out ", input.name))
    data_clust = read.table(data.clust, header = F, stringsAsFactors = F)
    split_pop = split(data_clust, data_clust[,3])
    pdf("./IBD_pop.pdf")
    par(mfrow = c(4,2))
    dev_pop = c()
    safe_pop = c()
    for(i in 1:length(split_pop)){
        write.table(split_pop[[i]][,c(1,2)], file = paste0(names(split_pop)[i], ".txt"), quote = F, row.names = F, col.names = F)
        system(paste0("plink --bfile ", input.name, " --keep ", names(split_pop)[i], ".txt --make-bed --out ", names(split_pop)[i]))
        system(paste0("plink --bfile ", names(split_pop)[i], " --extract ", input.name,".prune.in --genome --out ", names(split_pop)[i]))
        ibd_data = read.table(paste0(names(split_pop)[i], ".genome"), header = T)
        if(quantile(ibd_data$PI_HAT, 0.95) >= cutoff){
            dev_pop = rbind(dev_pop, c(names(split_pop)[i], round(mean(ibd_data$PI_HAT), 4), 
                                       round(median(ibd_data$PI_HAT), 4), round(quantile(ibd_data$PI_HAT, 0.95)[[1]], 4)))
        }else if(quantile(ibd_data$PI_HAT, 0.95) < cutoff){
            safe_pop = rbind(safe_pop, c(names(split_pop)[i], round(mean(ibd_data$PI_HAT), 4), 
                                        round(median(ibd_data$PI_HAT), 4), round(quantile(ibd_data$PI_HAT, 0.95)[[1]], 4)))
        }
        hist_info = hist(ibd_data$PI_HAT, plot = F)
        hist_info$prob = (hist_info$counts)/sum(hist_info$counts)
        plot(hist_info$mids, hist_info$prob, xlim = c(0,1), ylim = c(0,1), type = "h", col = "blue", xlab = "IBD", ylab = "Proportion of IBD", 
             main = paste0(names(split_pop)[i], " n=", length(unique(ibd_data[,1]))+1))
        abline(v = cutoff, lty = 2, col = "orange")
        abline(v = 0.5, lty = 2, col = "blue")
        abline(v = 0.25, lty = 2, col = "red")
    }
    dev.off()
    if(length(dev_pop) > 0){
        colnames(dev_pop) = c("POP", "MEAN", "MEDIAN", "95%")
        write.table(as.data.frame(dev_pop), file = "./dev_pop.txt", quote = F, row.names = F, col.names = T)
    }else{
        rm(dev_off)
    }
    if(length(safe_pop) > 0){
        colnames(safe_pop) = c("POP", "MEAN", "MEDIAN", "95%")
        write.table(as.data.frame(safe_pop), file = "./safe_pop.txt", quote = F, row.names = F, col.names = T)
    }else{
        rm(safe_pop)
    }
}

F_pop = function(input.name, data.clust){
    system(paste0("plink --bfile ",input.name," --indep-pairwise 50 5 0.2 --out ", input.name))
    system(paste0("plink --bfile ",input.name," --extract ",input.name,".prune.in --genome --out ", input.name))
    data_clust = read.table(data.clust, header = F, stringsAsFactors = F)
    split_pop = split(data_clust, data_clust[,3])
    pdf("./F.pdf")
    par(mfrow = c(4,2))
    for(i in 1:length(split_pop)){
        write.table(split_pop[[i]][,c(1,2)], file = paste0(names(split_pop)[i], ".txt"), quote = F, row.names = F, col.names = F)
        system(paste0("plink --bfile ", input.name, " --extract ", input.name,".prune.in --keep ", names(split_pop)[i], ".txt --make-bed --out ", names(split_pop)[i]))
        system(paste0("plink --bfile ", names(split_pop)[i], " --het --out ", names(split_pop)[i]))
        het_file = read.table(paste0(names(split_pop)[i], ".het"), header = T)
        b = hist(het_file$F, plot = F)
        b$prop = b$counts/nrow(het_file)
        plot(b$mids, b$prop, type = "h", xlim = c(-1,1), ylim = c(0,1),
             col = "blue", xlab = "F", ylab = "Proportion of F", 
             main = paste0(names(split_pop)[i], " n=", nrow(split_pop[[i]])))
    }
    dev.off()
}

prop_hom = function(input.name, data.clust){
    data_clust = read.table(data.clust, header = F, stringsAsFactors = F)
    split_pop = split(data_clust, data_clust[,3])
    pdf("./prop_hom.pdf")
    par(mfrow = c(4,2))
    for(i in 1:length(split_pop)){
        write.table(split_pop[[i]][,c(1,2)], file = paste0(names(split_pop)[i], ".txt"), quote = F, row.names = F, col.names = F)
        system(paste0("plink --bfile ", input.name, " --keep ", names(split_pop)[i], ".txt --make-bed --out ", names(split_pop)[i]))
        system(paste0("plink --bfile ", names(split_pop)[i], " --freqx --out ", names(split_pop)[i]))
        a = read.table(paste0(names(split_pop)[i], ".frqx"), header = T, sep = "\t")
        new_a = a %>% 
            filter(CHR %in% 1:22) %>%
            mutate(HOM.RATE = round(1-(C.HET./(nrow(split_pop[[i]])-C.MISSING.)), 3))
        b = hist(new_a$HOM.RATE, plot = F)
        b$prop = b$counts/nrow(new_a)
        plot(b$mids, b$prop, type = "h", xlim = c(0,1), ylim = c(0,0.5),
             col = "blue", xlab = "Homozygosity rate", ylab = "Proportion of HOM.RATE", 
             main = paste0(names(split_pop)[i], " n=", nrow(split_pop[[i]])))
        abline(v = 0.5, lty = 2, col = "red")
    }
    dev.off()
}

ind_qc_rm = function(input.name){
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
    fail_qc_inds = write.table(unique(a[order(a$V1),]), paste0(input.name, "_fail_ind_qc.txt"), row.names = F, col.names = F, quote = F)   
}

ind_qc_info = function(input.name){
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
        write.csv(df, file = paste0(input.name, "_ind_qc_info.csv"), quote = F, row.names = F)
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
        # pdf(paste0(input.name, "_threeVen.pdf"))
        # plot(threeVen)
        # dev.off()
    }else if(length(problem_all) == 4){
        AB = Reduce(intersect, c(problem_all[1], problem_all[3]))
        A = sub("\\..*", "", names(problem_all[1]))
        B = sub("\\..*", "", names(problem_all[3]))
        df = data.frame(AB)
        colnames(df) = paste(A, B, sep = "+")
        write.csv(df, file = paste0(input.name, "_ind_qc_info.csv"), quote = F, row.names = F)
        # twoVen = venneuler(c(A = length(problem_all[[1]]), B = length(problem_all[[3]]), "A&B" = length(out_table[[2]])))
        # twoVen$labels = c(
        #     paste0(A, "\n", length(problem_all[[1]])),
        #     paste0(B, "\n", length(problem_all[[3]]))
        # )
        # pdf(paste0(input.name, "_tewVen.pdf"))
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

missing_snp = function(input.name, missing_rate){
    system(paste0("plink --bfile ",input.name," --missing --out ", input.name))
    # Plot snp missing distribution
    # Load SNP frequency file and generate histogram
    b.frq = read.table(paste0(input.name, ".lmiss"), header = T)
    pdf("snpmiss_plot.pdf")
    plot(ecdf(b.frq$F_MISS), xlim = c(0,0.10), ylim = c(0,1), pch = 20, 
         main = "SNP Missingness Distribution", xlab = "Missingness Frequency", 
         ylab = "Fraction of SNPs", col = "blue")
    abline(v=missing_rate, col = "red", lty = 2)
    dev.off()
}

hwe_test = function(input.name, pvalue){
    system(paste0("plink --bfile ",input.name," --hardy --out ", input.name))
    hwe.file = read.table(paste0(input.name, ".hwe"), header = T)
    pdf("hwe_p_value.pdf")
    plot(ecdf(na.omit(hwe.file$P)), xlim = range(hwe.file$P,na.rm = T), 
         ylim = c(0,1), pch = 20, main = "Distribution of HWE P-Value", 
         xlab = "P-value", ylab = "Fraction of SNPs",col="blue")
    abline(v = pvalue, col = "red", lty = 2)
    dev.off()
}
