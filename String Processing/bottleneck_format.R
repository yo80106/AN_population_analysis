if(sum(installed.packages(.libPaths())[,1] %in% "stringr") == 0){
    install.packages("stringr")
    library(stringr)
}else{
    library(stringr)
}
args = commandArgs(trailingOnly = TRUE)

test = read.csv(args[1], header = FALSE, stringsAsFactors = FALSE)
head(test)
for(i in 2:ncol(test)){
    test[,i] = gsub("\\?", 00, test[,i])
    test[,i] = str_pad(test[,i], width =2, side="left", pad=0)
}
new = data.frame(test[,1])
for(t in seq(2,21,2)){
    new = cbind(new, paste0(test[,t],test[,t+1]))
}
colnames(new) = c("POP", 1:10)
new[,1] = paste0(test[,1],",")
write.table(new, file = args[2], sep = " ", quote = FALSE, row.names = FALSE)