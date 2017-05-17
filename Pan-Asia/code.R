setwd("~/Test/")
pan_fam = read.table("./Genotypes_All.fam", stringsAsFactors = F)
# Fill in Pan-Asia population list
library(stringr)
for(i in 1:nrow(pan_fam)){
    pop =  str_sub(pan_fam$V2[i],start = 1,end = 5)
    pan_fam$pop[i] = pop 
}
# Revice several population names
for(p in 1:nrow(pan_fam)){
    if(pan_fam$pop[p] == "CEU-N" | 
       pan_fam$pop[p] == "CHB-N" |
       pan_fam$pop[p] == "JPT-N" |
       pan_fam$pop[p] == "YRI-N"){
    pan_fam$pop[p] = substr(pan_fam$pop[p],1,3)
    }
}
write.table(pan_fam[,c(1,2,7)], file = "./data.clust", quote = FALSE, row.names = FALSE, col.names = FALSE)

system("plink --bfile Genotypes_All --freq --missing --within data.clust --out Genotypes_All")
system("gzip -c Genotypes_All.frq.strat > Genotypes_All.frq.strat.gz")
system("python plink2treemix.py Genotypes_All.frq.strat.gz Genotypes_All.treemix.gz")

## Bootstrap
setwd("~/Test/bootstrap/")
core = 5 # h is the number of cores that will be used
run = 20 # number of runs for each core; total number of bootstrap replicates is core * run  (20*5=100)
for (h in 1 : core){
    file.temp = paste("treemix_bootstrap_",h, ".sh", sep="")
    line.0 = paste("cd /home/chunyo/Test/bootstrap")
    cat(line.0, file= file.temp , append=TRUE,sep="\n")
    for (i in 1 : run){
        bootstrap.num <- (h-1)*run + i
        file.head=paste("shell_test", bootstrap.num, sep="-")
        line.1 = paste0("treemix -i Genotypes_All.treemix.gz -root YRI -bootstrap -k 75 -o Genotypes_All_", file.head)
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

sub_pop = c("ID-ML", "ID-DY", "PI-IR", "ID-TR", "PI-MA", "PI-UB", "ID-MT", "PI-UN", "PI-UI", "AX-AT", "AX-AM")
t_pop = lat.long[lat.long$Pop_name %in% sub_pop,]
t_pop = t_pop[-c(1,8),]
pdf(file = "./southeast_map.pdf")
map('world', xlim=c(90,160), ylim=c(-20,30))
points(t_pop$long.decimal, t_pop$lat.decimal, col = "red", pch = 20, cex = 1)
text(t_pop$long.decimal+4, t_pop$lat.decimal, labels = t_pop$Pop_name, cex = 0.7, col = "blue", font = 2)
dev.off()

lat.long = rbind(lat.long, list("Taitung", NA, NA, 22.45, 121.1))

pdf("./three_region_map.pdf")
austro = c("ID-ML", "ID-DY", "PI-IR", "ID-TR", "PI-MA", "PI-UB", "ID-MT", "PI-UN", "PI-UI")
tw = c("AX-AM", "AX-AT", "Taitung")
asiatic = c("ID-SU", "ID-JV", "MY-BD", "ID-JA", "MY-TM", "ID-TR", "ID-KR", "MY-MN", "ID-SB", "PI-AG")
mela = c("PI-MW", "PI-AT", "PI-AE", "ID-RA"," ID-SO", "ID-LA","ID-LE","ID-AL", "AX-ME")
austro_pop = lat.long[lat.long$Pop_name %in% austro,]
asiatic_pop = lat.long[lat.long$Pop_name %in% asiatic,]
mela_pop = lat.long[lat.long$Pop_name %in% mela,]
map('world', xlim=c(90,160), ylim=c(-20,30) ,fill = T, col = "gray")

austro_pop = austro_pop[-c(1,8),]
points(austro_pop$long.decimal, austro_pop$lat.decimal, col = "blue", pch = 20, cex = 1.3)
text(austro_pop$long.decimal+3, austro_pop$lat.decimal, labels = austro_pop$Pop_name, cex = 0.7, col = "blue", font = 2)

asiatic_pop = asiatic_pop[-c(9,12,13),]
points(asiatic_pop$long.decimal, asiatic_pop$lat.decimal, col = "dark green", pch = 20, cex = 1.3)
text(asiatic_pop$long.decimal+3, asiatic_pop$lat.decimal, labels = asiatic_pop$Pop_name, cex = 0.7, col = "dark green", font = 2)

mela_pop = mela_pop[-c(7:10,12:18,20:23),]
points(mela_pop$long.decimal, mela_pop$lat.decimal, col = "orange", pch = 20, cex = 1.3)
text(mela_pop$long.decimal+3, mela_pop$lat.decimal, labels = mela_pop$Pop_name, cex = 0.7, col = "orange", font = 2)
dev.off()

pdf("./austro_map.pdf")
austro_pop = lat.long[lat.long$Pop_name %in% austro,]
map('world', xlim=c(90,160), ylim=c(-20,30) ,fill = T, col = "gray")
austro_pop = austro_pop[-c(1,8),]
points(austro_pop$long.decimal, austro_pop$lat.decimal, col = "blue", pch = 20, cex = 1.3)
text(austro_pop$long.decimal+3, austro_pop$lat.decimal, labels = austro_pop$Pop_name, cex = 0.7, col = "blue", font = 2)
dev.off()

pdf("./austro_map_diff_col.pdf")
austro_pop = lat.long[lat.long$Pop_name %in% austro,]
map('world', xlim=c(90,160), ylim=c(-20,30) ,fill = T, col = "gray")
austro_pop = austro_pop[-6,]
points(austro_pop$long.decimal, austro_pop$lat.decimal, col = "orange", pch = 20, cex = 1.3)
text(austro_pop$long.decimal+3, austro_pop$lat.decimal, labels = austro_pop$Pop_name, cex = 0.7, col = "orange", font = 2)
tw_pop = lat.long[lat.long$Pop_name %in% tw,]
points(tw_pop$long.decimal, tw_pop$lat.decimal, col = "blue", pch = 20, cex = 1.3)
text(tw_pop$long.decimal+3, tw_pop$lat.decimal, labels = tw_pop$Pop_name, cex = 0.7, col = "blue", font = 2)
dev.off()