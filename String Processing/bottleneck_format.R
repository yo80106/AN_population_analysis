test = read.csv("~/Downloads/cw22type123.csv", header = FALSE, stringsAsFactors = FALSE)
str(test)
head(test)
library(stringr)
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
write.table(new, file = "~/Desktop/angela.txt", sep = " ", quote = FALSE, row.names = FALSE)