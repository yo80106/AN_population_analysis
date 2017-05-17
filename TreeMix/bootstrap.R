setwd("~/TreeMix/test_bootstrap/")

# Load sample data
clust = read.table("./data.clust")
test_fam = read.table("./austronesian.fam")
clust = clust[clust$V2 %in% test_fam$V2,]
write.table(clust,"./austro.clust",row.names = F,col.names = F,quote = F)
    
system("plink --bfile austronesian --freq --missing --within austro.clust --out austro_bootstrap")
system("gzip -c austro_bootstrap.frq.strat > austro_bootstrap.frq.strat.gz")
system("python plink2treemix.py austro_bootstrap.frq.strat.gz austro_bootstrap.treemix.gz")

for (n in 1:50) {
    system(paste("treemix -i test_bootstrap.treemix.gz -bootstrap -k 1000 -o test",n,sep = "")) 
}

system("mkdir bootstrap_data")
system("mv *.treeout.gz bootstrap_data/")

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

# Mutiple Cores
a = c()
for (n in 1:50) {
   a[n] = paste("treemix -i austro_bootstrap.treemix.gz -root San -bootstrap -k 1000 -o austro",n,sep = "") 
}
lapply(a, system)
library(parallel)
?mclapply

# Load real data 
library(parallel)
a = c()
for (n in 1:20) {
    a[n] = paste0("taskset -cp ", n," treemix -i austro_bootstrap.treemix.gz -bootstrap -k 500 -o test",n) 
}
cat(a, file = "./command.sh", sep = "\n")
lapply(a, system)
mclapply(a, system, mc.cores = 2)
cl = makeCluster(10)
parLapply(cl, a, system)

