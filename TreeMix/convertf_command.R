setwd("~/AdmixTools/data/")
# The syntax of convertf is "../bin/convertf -p parfile".
# genotypename:    example.ped
# snpname:         example.pedsnp # or example.map, either works 
# indivname:       example.pedind # or example.ped, either works
# outputformat:    EIGENSTRAT
# genotypeoutname: example.eigenstratgeno
# snpoutname:      example.snp
# indivoutname:    example.ind
# familynames:     NO

system("mv ./data/final_merge_no_LD.bim ./data/final_merge_no_LD.pedsnp")
system("mv ./data/final_merge_no_LD.fam ./data/final_merge_no_LD.pedind")
parfile = c("genotypename:  final_merge_no_LD.ped",
            "snpname:   final_merge_no_LD.pedsnp",
            "indivname: final_merge_no_LD.pedind",
            "outputformat:  EIGENSTRAT",
            "genotypeoutname:   final_merge_no_LD.geno",
            "snpoutname:    final_merge_no_LD.snp",
            "indivoutname:  final_merge_no_LD.ind",
            "familynames:   NO")
cat(parfile, file = "./parfile", sep = "\n")
system(" ~/AdmixTools/bin/convertf -p parfile")

pop = read.table("./data.clust")
ind = read.table("./final_merge_no_LD.ind")
sum(!pop$V2 %in% ind$V1)
ind$V3 = pop$V3
head(ind)
write.table(ind, file = "./new_final_merge_no_LD.ind",quote = F,col.names = F,row.names = F)

## IN SERVER
# MixMapper SECTION
system("../compute_moment_stats new_final_merge_no_LD.ind final_merge_no_LD.snp final_merge_no_LD.geno test n 50 15 > compute_moment_stats_stdout.txt 2> compute_moment_stats_stderr.txt")

# MixMapper SECTION selected populations
system("../compute_moment_stats new_final_merge_no_LD.ind final_merge_no_LD.snp final_merge_no_LD.geno test n 100 15 > compute_moment_stats_stdout.txt 2> compute_moment_stats_stderr.txt")