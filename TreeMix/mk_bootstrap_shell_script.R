core = 25 # h is the number of cores that will be used
run = 4 # number of runs for each core; total number of bootstrap replicates is core * run  (20*5=100)
for (h in 1 : core){
	file.temp = paste("treemix_bootstrap_",h, ".sh", sep="")
	line.0 = paste("cd /home/chunyo/TreeMix/bootstrap/bootstrap_data")
	cat(line.0, file= file.temp , append=TRUE,sep="\n")
	for (i in 1 : run){
		bootstrap.num <- (h-1)*run + i
		file.head=paste("shell_test", bootstrap.num, sep="-")
		line.1 = paste("treemix -i final_merge_no_LD.treemix.gz -root San -bootstrap -k 15 -o", file.head, sep=" ")
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