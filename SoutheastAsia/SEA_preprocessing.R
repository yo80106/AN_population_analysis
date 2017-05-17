setwd("~/SoutheastAsia/")
sea_fam = read.table("./SEA_730K.fam", header = FALSE)
dim(sea_fam)

############################# Taitung + SoutheastAsia ##########################
sea_bim = read.table("./SEA_730K.bim", header = FALSE, stringsAsFactors = FALSE)
str(sea_bim)
taitung_bim = read.table("./new_clean_per_snp.bim", header = FALSE, stringsAsFactors = FALSE)
str(taitung_bim)
t_merge = merge(taitung_bim, sea_bim, by.x = "V2", by.y = "V2")
dim(t_merge)
head(t_merge)
write.table(t_merge[,1], file = "./merge_snp_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile new_clean_per_snp --extract merge_snp_list.txt --make-bed --out taitung_extract")
system("plink --bfile SEA_730K --extract merge_snp_list.txt --make-bed --out southeastasia_extract")
system("plink --bfile taitung_extract --bmerge southeastasia_extract --make-bed --out test1")
system("plink --bfile taitung_extract --flip test1-merge.missnp --make-bed --out taitung_extract_flip")
system("plink --bfile taitung_extract_flip --bmerge southeastasia_extract --make-bed --out test2")
system("plink --bfile taitung_extract_flip --exclude test2-merge.missnp --make-bed --out exclude_taitung_ex_flip")
system("plink --bfile exclude_taitung_ex_flip --bmerge southeastasia_extract --make-bed --out test3")
system("mv test3.bed taitung_asia_merge.bed")
system("mv test3.bim taitung_asia_merge.bim")
system("mv test3.fam taitung_asia_merge.fam")

merge_fam = read.table("./taitung_asia_merge.fam", stringsAsFactors = FALSE)
taitung_list = read.csv("~/PCA/R_analysis/96samples_sex_list.csv",header = F, stringsAsFactors = FALSE)

for(i in 1:nrow(merge_fam)){
    pop =  gsub("\\d", "", merge_fam$V2[i])
    merge_fam$pop[i] = pop 
}

for(t in 1:nrow(merge_fam)){
    for(k in 1:nrow(taitung_list)){
        if(merge_fam$V2[t] == taitung_list$V2[k]){
            merge_fam$pop[t] = taitung_list$V4[k]
        }
    }
}

source("~/Git/PCA/pca_function.R")
data.clust1 = write.table(merge_fam[,c(1,2,7)], file = "./data.clust1", quote = FALSE, col.names = FALSE, row.names = FALSE)
pop.list1 = write.table(unique(merge_fam$pop), file = "./pop.list1", quote = FALSE, col.names = FALSE, row.names = FALSE)
my_pca("taitung_asia_merge", "~/SoutheastAsia/data.clust1", "~/SoutheastAsia/pop.list1", "taitung_asia_merge")

merge_fam[merge_fam$V2 == "TDC480",]
merge_fam[merge_fam$V2 == "TDC480_2",]
merge_fam = merge_fam[-c(88,95),]
table(merge_fam$pop)
data.clust2 = write.table(merge_fam[,c(1,2,7)], file = "./data.clust2", quote = FALSE, col.names = FALSE, row.names = FALSE)
my_pca("taitung_asia_merge", "~/SoutheastAsia/data.clust2", "~/SoutheastAsia/pop.list1", "taitung_asia_merge_no_outlier")

############################ Taitung + HGDP + SoutheastAsia#####################
setwd("~/SoutheastAsia/Taitung_Asia_HGDP/")
hgdp_bim = read.table("./HGDP_Map.bim", stringsAsFactors = FALSE)
str(hgdp_bim)
taitung_asia_bim = read.table("./taitung_asia_merge.bim", stringsAsFactors = FALSE)
str(taitung_asia_bim)
t_merge = merge(taitung_asia_bim, hgdp_bim, by.x = "V2", by.y = "V2")
str(t_merge)
write.table(t_merge[,1], file = "./taitung_asia_hgdp_snp_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile taitung_asia_merge --extract taitung_asia_hgdp_snp_list.txt --make-bed --out taitung_asia_extract")
system("plink --bfile HGDP_Map --extract taitung_asia_hgdp_snp_list.txt --make-bed --out hgdp_extract")
system("plink --bfile taitung_asia_extract --bmerge hgdp_extract --make-bed --out test5")
system("plink --bfile taitung_asia_extract --flip test5-merge.missnp --make-bed --out taitung_asia_extract_flip")
system("plink --bfile taitung_asia_extract_flip --bmerge hgdp_extract --make-bed --out test6")
system("mv test6.bed taitung_asia_hgdp.bed")
system("mv test6.bim taitung_asia_hgdp.bim")
system("mv test6.fam taitung_asia_hgdp.fam")

taitung_asia_hgdp_fam = read.table("./taitung_asia_hgdp.fam", stringsAsFactors = FALSE)
str(taitung_asia_hgdp_fam)
sea_list = read.csv("~/SoutheastAsia/SEA_ind_info.csv", header = FALSE, stringsAsFactors = FALSE)
str(sea_list)
final_sample_id = read.table("~/SoutheastAsia/final_data.clust", header = FALSE, stringsAsFactors = FALSE)
str(final_sample_id)

taitung_asia_hgdp_fam$pop = NA
taitung_asia_hgdp_fam$pop[!is.na(match(taitung_asia_hgdp_fam[,2], sea_list[,1]))] = sea_list[,2][na.omit(match(taitung_asia_hgdp_fam[,2], sea_list[,1]))]
taitung_asia_hgdp_fam$pop[!is.na(match(taitung_asia_hgdp_fam[,2], final_sample_id[,2]))] = final_sample_id[,3][na.omit(match(taitung_asia_hgdp_fam[,2], final_sample_id[,2]))]
write.table(taitung_asia_hgdp_fam[,c(1,2,7)], file = "./data.clust", quote = FALSE, row.names = FALSE, col.names = FALSE)

## TreeMix
system("plink --bfile taitung_asia_hgdp --freq --missing --within data.clust --out taitung_asia_hgdp")
system("gzip -c taitung_asia_hgdp.frq.strat > taitung_asia_hgdp.frq.strat.gz")
system("python plink2treemix.py taitung_asia_hgdp.frq.strat.gz taitung_asia_hgdp.treemix.gz")

n = "taitung_asia_hgdp"
file.temp = paste0("treemix_bootstrap_",n, ".sh")
line.0 = paste0("cd /home/chunyo/SoutheastAsia")
line.1 = paste0("treemix -i taitung_asia_hgdp.treemix.gz -root San -bootstrap -k 400 -o taitung_asia_hgdp")
cat(line.0, file= file.temp , append=TRUE,sep="\n")
cat(line.1, file= file.temp , append=TRUE,sep="\n")

file.temp = paste0("pbs_treemix_bootstrap_",n,".sh")
line.1 = paste0("#PBS -N treemix_bootstrap_",n)
line.2 = paste0("#PBS -e treemix_bootstrap_",n,".err")
line.3 = paste0("#PBS -o treemix_bootstrap_",n,".out")
line.4 = paste0("cd $PBS_O_WORKDIR")
line.5 = paste0("sh treemix_bootstrap_",n,".sh")
cat(line.1, file= file.temp , append=TRUE,sep="\n")
cat(line.2, file= file.temp , append=TRUE,sep="\n")
cat(line.3, file= file.temp , append=TRUE,sep="\n")
cat(line.4, file= file.temp , append=TRUE,sep="\n")
cat(line.5, file= file.temp , append=TRUE,sep="\n")

# qsub pbs_treemix_bootstrap_taitung_asia_hgdp.sh
## Add migration
n = "taitung_asia_hgdp"
for(i in 2:10){
    file.temp = paste0("treemix_bootstrap_",n,"_m", i, ".sh")
    line.0 = paste0("cd /home/chunyo/SoutheastAsia/Taitung_Asia_HGDP/migration")
    line.1 = paste0("treemix -i taitung_asia_hgdp.treemix.gz -root San -bootstrap -k 400 -m ", i," -o taitung_asia_hgdp_m", i)
    cat(line.0, file= file.temp , append=TRUE,sep="\n")
    cat(line.1, file= file.temp , append=TRUE,sep="\n")
}

for(t in 2:10){
    file.temp2 = paste0("pbs_treemix_bootstrap_",n, "_m", t, ".sh")
    line.1 = paste0("#PBS -N treemix_bootstrap_", n, "_m", t)
    line.2 = paste0("#PBS -e treemix_bootstrap_", n, "_m", t, ".err")
    line.3 = paste0("#PBS -o treemix_bootstrap_", n, "_m", t, ".out")
    line.4 = paste0("cd $PBS_O_WORKDIR")
    line.5 = paste0("sh treemix_bootstrap_",n, "_m", t, ".sh")
    cat(line.1, file= file.temp2 , append=TRUE,sep="\n")
    cat(line.2, file= file.temp2 , append=TRUE,sep="\n")
    cat(line.3, file= file.temp2 , append=TRUE,sep="\n")
    cat(line.4, file= file.temp2 , append=TRUE,sep="\n")
    cat(line.5, file= file.temp2 , append=TRUE,sep="\n")
}
source("~/TreeMix/plotting_funcs.R")
setwd("~/SoutheastAsia/Taitung_Asia_HGDP/TreeMix/migration/")
for(z in 2:10){
    pdf(paste0("./taitung_asia_hgdp_m", z, ".pdf"))
    plot_tree(paste0("taitung_asia_hgdp_m", z))
    dev.off()
}

## Bootstrap
setwd("~/SoutheastAsia/Taitung_Asia_HGDP/TreeMix/Bootstrap/")
core = 10 # h is the number of cores that will be used
run = 10 # number of runs for each core; total number of bootstrap replicates is core * run  (20*5=100)
for (h in 1 : core){
    file.temp = paste("treemix_bootstrap_",h, ".sh", sep="")
    line.0 = paste("cd /home/chunyo/SoutheastAsia/Taitung_Asia_HGDP/Bootstrap")
    cat(line.0, file= file.temp , append=TRUE,sep="\n")
    for (i in 1 : run){
        bootstrap.num <- (h-1)*run + i
        file.head=paste("shell_test", bootstrap.num, sep="-")
        line.1 = paste0("treemix -i taitung_asia_hgdp.treemix.gz -root San -bootstrap -k 400 -o taitung_asia_hgdp_", file.head)
        cat(line.1, file= file.temp , append=TRUE,sep="\n")
    }
}

for(n in 1:core){
    file.temp = paste0("pbs_treemix_bootstrap_",n,".sh")
    line.1 = paste0("#PBS -N treemix_bootstrap_",n)
    line.2 = paste0("#PBS -e treemix_bootstrap_",n,".err")
    line.3 = paste0("#PBS -o treemix_bootstrap_",n,".out")
    line.4 = paste0("cd $PBS_O_WORKDIR")
    line.5 = paste0("sh treemix_bootstrap_",n,".sh")
    cat(line.1, file= file.temp , append=TRUE,sep="\n")
    cat(line.2, file= file.temp , append=TRUE,sep="\n")
    cat(line.3, file= file.temp , append=TRUE,sep="\n")
    cat(line.4, file= file.temp , append=TRUE,sep="\n")
    cat(line.5, file= file.temp , append=TRUE,sep="\n")
}
setwd("~/SoutheastAsia/Taitung_Asia_HGDP/TreeMix/Bootstrap/treeout/")
file_list = list.files("./",pattern = ".gz")
for(x in file_list){
    if(!exists("read_line")){
        read_line = readLines(x)
        next()
    }
    if(exists("read_line")){
        temp = readLines(x)
        read_line = rbind(read_line,temp,deparse.level = 0)
        rm(temp)
    }
}
concat = gsub('[0-9:.]+', '', read_line)
concat = gsub("e-","",concat)
write.table(concat,"./newick",row.names = F,col.names = F,quote = F)

## Dstat
setwd("~/SoutheastAsia/Taitung_Asia_HGDP/Dstat/")
dstat.1 = data.frame(pop1 = rep(c("Amis", "Bunun", "Paiwan", "Puyuma"), each = 3),
                     pop2 = rep(c("Murut", "Filipino", "Dusun"), 4),
                     pop3 = rep("Kankanaey", 3*4),
                     pop4 = rep("Yoruba", 3*4))
write.table(dstat.1, file = "./dstat_pop_1", quote = FALSE, row.names = FALSE, col.names = FALSE)

system("plink --bfile taitung_asia_hgdp --recode --out taitung_asia_hgdp")
# The syntax of convertf is "../bin/convertf -p parfile".
# genotypename:    example.ped
# snpname:         example.pedsnp # or example.map, either works 
# indivname:       example.pedind # or example.ped, either works
# outputformat:    EIGENSTRAT
# genotypeoutname: example.eigenstratgeno
# snpoutname:      example.snp
# indivoutname:    example.ind
# familynames:     NO

system("mv ./taitung_asia_hgdp.bim ./taitung_asia_hgdp.pedsnp")
system("mv ./taitung_asia_hgdp.fam ./taitung_asia_hgdp.pedind")
parfile = c("genotypename:  taitung_asia_hgdp.ped",
            "snpname:   taitung_asia_hgdp.pedsnp",
            "indivname: taitung_asia_hgdp.pedind",
            "outputformat:  EIGENSTRAT",
            "genotypeoutname:   taitung_asia_hgdp.geno",
            "snpoutname:    taitung_asia_hgdp.snp",
            "indivoutname:  taitung_asia_hgdp.ind",
            "familynames:   NO")
cat(parfile, file = "./convert_parfile", sep = "\n")
system(" ~/AdmixTools/bin/convertf -p convert_parfile")

data_clust = read.table("~/SoutheastAsia/Taitung_Asia_HGDP/data.clust", stringsAsFactors = FALSE)
str(data_clust)
ind_file = read.table("./taitung_asia_hgdp.ind", stringsAsFactors = FALSE)
str(ind_file)
ind_file$V3[!is.na(match(ind_file[,1], data_clust[,2]))] = data_clust[,3][na.omit(match(ind_file[,1], data_clust[,2]))]
write.table(ind_file, file = "./taitung_asia_hgdp.ind", quote = FALSE, row.names = FALSE, col.names = FALSE)

parfile = c("genotypename:  taitung_asia_hgdp.geno",
            "snpname:   taitung_asia_hgdp.snp",
            "indivname: taitung_asia_hgdp.ind",
            "popfilename: dstat_pop_1")
cat(parfile, file = "./dstat_par1", sep = "\n")
system("~/AdmixTools/bin/qpDstat -p dstat_par1 > dstat_par1_out")

################# Taitung + HGDP + Pan-Asia + SoutheastAsia ####################
final_bim = read.table("./final_merge_no_LD.bim", stringsAsFactors = FALSE)
str(final_bim)
sea_bim = read.table("./SEA_730K.bim", header = FALSE, stringsAsFactors = FALSE)
str(sea_bim)

t_merge = merge(final_bim, sea_bim, by.x = "V2", by.y = "V2")
dim(t_merge)
head(t_merge)
write.table(t_merge[,1], file = "./final_sea_snp_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

system("plink --bfile final_merge_no_LD --extract final_sea_snp_list.txt --make-bed --out final_extract")
system("plink --bfile SEA_730K --extract final_sea_snp_list.txt --make-bed --out southeastasia_extract2")
system("plink --bfile final_extract --bmerge southeastasia_extract2 --make-bed --out test4")
system("plink --bfile final_extract --flip test4-merge.missnp --make-bed --out final_extract_flip")
system("plink --bfile final_extract_flip --bmerge southeastasia_extract2 --make-bed --out test5")
system("mv test5.bed final_asia_merge.bed")
system("mv test5.bim final_asia_merge.bim")
system("mv test5.fam final_asia_merge.fam")

final_asia_fam = read.table("~/SoutheastAsia/final_asia_merge.fam", stringsAsFactors = FALSE)
str(final_asia_fam)
final_sample_id = read.table("~/SoutheastAsia/final_data.clust", header = FALSE, stringsAsFactors = FALSE)
str(final_sample_id)
sea_list = read.csv("~/SoutheastAsia/SEA_ind_info.csv", header = FALSE, stringsAsFactors = FALSE)
str(sea_list)

final_asia_fam$pop = NA
final_asia_fam$pop[!is.na(match(final_asia_fam[,2], final_sample_id[,2]))] = final_sample_id[,3][na.omit(match(final_asia_fam[,2], final_sample_id[,2]))]
final_asia_fam$pop[!is.na(match(final_asia_fam[,2], sea_list[,1]))] = sea_list[,2][na.omit(match(final_asia_fam[,2], sea_list[,1]))]
write.table(final_asia_fam[,c(1,2,7)], file = "./data.clust3", quote = FALSE, row.names = FALSE, col.names = FALSE)

## TreeMix
setwd("~/SoutheastAsia/")
system("plink --bfile final_asia_merge --freq --missing --within data.clust3 --out final_asia_merge")
system("gzip -c final_asia_merge.frq.strat > final_asia_merge.frq.strat.gz")
system("python plink2treemix.py final_asia_merge.frq.strat.gz final_asia_merge.treemix.gz")

n = "final_asia"
file.temp = paste0("treemix_bootstrap_",n, ".sh")
line.0 = paste0("cd /home/chunyo/SoutheastAsia")
line.1 = paste0("treemix -i final_asia_merge.treemix.gz -root San -bootstrap -k 15 -o final_asia_merge")
cat(line.0, file= file.temp , append=TRUE,sep="\n")
cat(line.1, file= file.temp , append=TRUE,sep="\n")

file.temp = paste0("pbs_treemix_bootstrap_",n,".sh")
line.1 = paste0("#PBS -N treemix_bootstrap_",n)
line.2 = paste0("#PBS -e treemix_bootstrap_",n,".err")
line.3 = paste0("#PBS -o treemix_bootstrap_",n,".out")
line.4 = paste0("cd $PBS_O_WORKDIR")
line.5 = paste0("sh treemix_bootstrap_",n,".sh")
cat(line.1, file= file.temp , append=TRUE,sep="\n")
cat(line.2, file= file.temp , append=TRUE,sep="\n")
cat(line.3, file= file.temp , append=TRUE,sep="\n")
cat(line.4, file= file.temp , append=TRUE,sep="\n")
cat(line.5, file= file.temp , append=TRUE,sep="\n")

# qsub pbs_treemix_bootstrap_final_asia.sh
## Bootstrap
setwd("~/SoutheastAsia/TreeMix/Bootstrap/")
core = 25 # h is the number of cores that will be used
run = 4 # number of runs for each core; total number of bootstrap replicates is core * run  (20*5=100)
for (h in 1 : core){
    file.temp = paste("treemix_bootstrap_",h, ".sh", sep="")
    line.0 = paste("cd /home/chunyo/SoutheastAsia/Bootstrap")
    cat(line.0, file= file.temp , append=TRUE,sep="\n")
    for (i in 1 : run){
        bootstrap.num <- (h-1)*run + i
        file.head=paste("shell_test", bootstrap.num, sep="-")
        line.1 = paste0("treemix -i final_asia_merge.treemix.gz -root San -bootstrap -k 15 -o final_asia_merge_", file.head)
        cat(line.1, file= file.temp , append=TRUE,sep="\n")
    }
}

for(n in 1:core){
    file.temp = paste0("pbs_treemix_bootstrap_",n,".sh")
    line.1 = paste0("#PBS -N treemix_bootstrap_",n)
    line.2 = paste0("#PBS -e treemix_bootstrap_",n,".err")
    line.3 = paste0("#PBS -o treemix_bootstrap_",n,".out")
    line.4 = paste0("cd $PBS_O_WORKDIR")
    line.5 = paste0("sh treemix_bootstrap_",n,".sh")
    cat(line.1, file= file.temp , append=TRUE,sep="\n")
    cat(line.2, file= file.temp , append=TRUE,sep="\n")
    cat(line.3, file= file.temp , append=TRUE,sep="\n")
    cat(line.4, file= file.temp , append=TRUE,sep="\n")
    cat(line.5, file= file.temp , append=TRUE,sep="\n")
}

setwd("~/SoutheastAsia/Taitung_Asia_HGDP/TreeMix/Bootstrap/treeout/")
file_list = list.files("./",pattern = ".gz")
for(x in file_list){
    if(!exists("read_line")){
        read_line = readLines(x)
        next()
    }
    if(exists("read_line")){
        temp = readLines(x)
        read_line = rbind(read_line,temp,deparse.level = 0)
        rm(temp)
    }
}
concat = gsub('[0-9:.]+', '', read_line)
concat = gsub("e-","",concat)
write.table(concat,"./newick",row.names = F,col.names = F,quote = F)


## Pop3 test
setwd("~/SoutheastAsia/Pop3Test/")
system("plink --bfile final_asia_merge --recode --out final_asia_merge")
# The syntax of convertf is "../bin/convertf -p parfile".
# genotypename:    example.ped
# snpname:         example.pedsnp # or example.map, either works 
# indivname:       example.pedind # or example.ped, either works
# outputformat:    EIGENSTRAT
# genotypeoutname: example.eigenstratgeno
# snpoutname:      example.snp
# indivoutname:    example.ind
# familynames:     NO

system("mv ./final_asia_merge.bim ./final_asia_merge.pedsnp")
system("mv ./final_asia_merge.fam ./final_asia_merge.pedind")
parfile = c("genotypename:  final_asia_merge.ped",
            "snpname:   final_asia_merge.pedsnp",
            "indivname: final_asia_merge.pedind",
            "outputformat:  EIGENSTRAT",
            "genotypeoutname:   final_asia_merge.geno",
            "snpoutname:    final_asia_merge.snp",
            "indivoutname:  final_asia_merge.ind",
            "familynames:   NO")
cat(parfile, file = "./parfile", sep = "\n")
system(" ~/AdmixTools/bin/convertf -p parfile")

parfile = read.table("./final_asia_merge.ind", stringsAsFactors = FALSE)
str(parfile)
data_clust = read.table("~/SoutheastAsia/data.clust3", stringsAsFactors = FALSE)
str(data_clust)
parfile$V3[!is.na(match(parfile[,1], data_clust[,2]))] = data_clust[,3][na.omit(match(parfile[,1], data_clust[,2]))]
write.table(parfile, file = "./final_asia_merge.ind", quote = FALSE, row.names = FALSE, col.names = FALSE)

pop3.1 = data.frame(tested = rep(c("Amis", "AX-AM", "AX-AT", "Bunun", "Paiwan", "Puyuma"), each = 8),
                    ref1 = rep("CHB", 6*8),
                    ref2 = rep(c("Kankanaey", "PI-UB", "PI-MA", "Filipino", "ID-MT", "ID-TR", "PI-UN", "PI-UI"), 6))
write.table(pop3.1, file = "./pop3test_pop_1", quote = FALSE, row.names = FALSE, col.names = FALSE)

# DESCRIPTION OF EACH PARAMETER in parfile:
#     
# genotypename:   input genotype file (in eigenstrat format)
# snpname:   input snp file      (in eigenstrat format)
# indivname:   input indiv file    (in eigenstrat format)
# popfilename:  list_qp3test (contains 3 populations on each line <Source1 (A)> <Source2 (B)> < Target (C)>

parfile = c("genotypename:  final_asia_merge.geno",
            "snpname:   final_asia_merge.snp",
            "indivname: final_asia_merge.ind",
            "popfilename: pop3test_pop_1")
cat(parfile, file = "./qp3pop_par1", sep = "\n")
system("~/AdmixTools/bin/qp3Pop -p qp3pop_par1 >qp3pop_par1_out")

## DStats
setwd("~/SoutheastAsia/Dstats/")
dstat.1 = data.frame(pop1 = rep(c("Amis", "AX-AM", "AX-AT", "Bunun", "Paiwan", "Puyuma"), each = 8),
                    pop2 = rep("CHB", 6*8),
                    pop3 = rep(c("Kankanaey", "PI-UB", "PI-MA", "Filipino", "ID-MT", "ID-TR", "PI-UN", "PI-UI"), 6),
                    pop4 = rep("Yoruba", 6*8))
write.table(dstat.1, file = "./dstat_pop_1", quote = FALSE, row.names = FALSE, col.names = FALSE)
parfile = c("genotypename:  final_asia_merge.geno",
            "snpname:   final_asia_merge.snp",
            "indivname: final_asia_merge.ind",
            "popfilename: dstat_pop_1")
cat(parfile, file = "./dstat_par1", sep = "\n")
system("~/AdmixTools/bin/qpDstat -p dstat_par1 > dstat_par1_out")

