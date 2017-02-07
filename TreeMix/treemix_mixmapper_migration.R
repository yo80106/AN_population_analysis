setwd("~/TreeMix/select_pop_migration/")
## TreeMix SECTION
# Select populations
data.clust = read.table("./all.data.clust", stringsAsFactors = FALSE)
new.clust = data.clust[data.clust$V3 %in% c("Paiwan","Amis","Bunun","Puyuma","AX-AM","AX-AT",
                                "Miaozu","She","CN-JI","Lahu","CN-WA","Yizu","Naxi",
                                "TH-HM","TH-PP","TH-TN","TH-PL","Karitiana","Surui",
                                "Papuan","Mandenka","Yoruba","PI-AG","PI-AT","PI-AE",
                                "PI-IR","PI-MW","PI-MA","PI-UN","PI-UI","ID-AL","ID-SB",
                                "ID-LA","ID-LE","ID-SO","ID-RA","ID-TR","ID-TB","MY-BD",
                                "ID-DY","ID-JA","ID-JV","ID-ML","SG-MY","ID-SU","ID-KR",
                                "MY-KN","MY-MN","ID-MT","PI-UB","MY-TM","NANMelanesian",
                                "CN-CC","CN-JN","CN-GA","CN-HM","TH-TL","TH-TU"),]
write.table(new.clust, file = "./select.data.clust",row.names = FALSE,col.names = FALSE,quote = FALSE)
system("plink --bfile final_merge_no_LD --freq --missing --within select.data.clust --out select_pop")
system("gzip -c select_pop.frq.strat > select_pop.frq.strat.gz")
system("python plink2treemix.py select_pop.frq.strat.gz select_pop.treemix.gz")

# Results
setwd("~/TreeMix/select_pop_migration/TreeMix_results/")
source("~/TreeMix/plotting_funcs.R")
par(mfrow = c(4,2), mar = c(3.5, 3.5, 2, 1))
pdf("pasr1.pdf")
for(i in 2:20){
    plot_tree(paste0("mig_",i))
}
dev.off()

# IN SERVER
mig = 20
for (h in 2 : mig){
    file.temp = paste0("mig_",h,".sh")
    line.0 = paste0("cd /home/chunyo/TreeMix/migration")
    line.1 = paste0("treemix -i select_pop.treemix.gz -root Yoruba -m ",h," -o mig_",h)
    cat(line.0, file= file.temp , append=TRUE,sep="\n")
    cat(line.1, file= file.temp , append=TRUE,sep="\n")
}

for(n in 2:mig){
    file.temp = paste0("pbs_mig_",n,".sh")
    line.1 = paste0("#PBS -N mig_",n)
    line.2 = paste0("#PBS -e mig_",n,".err")
    line.3 = paste0("#PBS -o mig_",n,".out")
    line.4 = paste0("cd $PBS_O_WORKDIR")
    line.5 = paste0("sh mig_",n,".sh")
    cat(line.1, file= file.temp , append=TRUE,sep="\n")
    cat(line.2, file= file.temp , append=TRUE,sep="\n")
    cat(line.3, file= file.temp , append=TRUE,sep="\n")
    cat(line.4, file= file.temp , append=TRUE,sep="\n")
    cat(line.5, file= file.temp , append=TRUE,sep="\n")
}
## MixMapper SECTION
setwd("./Mixmapper/")
ind = data.frame(new.clust$V1, new.clust$V2)
write.table(ind, file = "ind.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("plink --bfile final_merge_no_LD --keep ind.txt --make-bed --out select_pop")
system("plink --bfile select_pop --recode --out select_pop")
system("mv select_pop.bim select_pop.pedsnp")
system("mv select_pop.fam select_pop.pedind")

parfile = c("genotypename:  select_pop.ped",
            "snpname:   select_pop.pedsnp",
            "indivname: select_pop.pedind",
            "outputformat:  EIGENSTRAT",
            "genotypeoutname:   select_pop.geno",
            "snpoutname:    select_pop.snp",
            "indivoutname:  select_pop.ind",
            "familynames:   NO")
cat(parfile, file = "./parfile", sep = "\n")
system(" ~/AdmixTools/bin/convertf -p parfile")

ind = read.table("select_pop.ind")
sum(!new.clust$V2 %in% ind$V1)
ind$V3 = new.clust$V3
head(ind)
write.table(ind, file = "./new_select_pop.ind",quote = F,col.names = F,row.names = F)

# IN SERVER
system("/home/chunyo/MixMapper/MixMapper_v2.0/compute_moment_stats new_select_pop.ind select_pop.snp select_pop.geno mig n 100 15 > compute_moment_stats_stdout.txt 2> compute_moment_stats_stderr.txt")
setwd("./results/")
stdout = readLines("compute_moment_stats_stdout.txt")
new.stdout = as.data.frame(stdout[!grepl("number of", stdout)])
write.table(new.stdout, file = "new_stdout.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("/home/chunyo/MixMapper/MixMapper_v2.0/compute_most_additive_trees mig.f2.tab 10000 new_stdout.txt > most_additive_trees.txt")