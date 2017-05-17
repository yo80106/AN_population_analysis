setwd("~/Desktop/")
inp_list = list.files("./", pattern = ".inp")
for(n in 1:22){
  sour = readLines("./NCAM_flipped_MAF0.05_chr",n,"_rehh.inp")
}
pas = readLines("./FK.txt")
pas.sep = strsplit(pas, "  ")
pas.vec = unlist(lapply(pas.sep, function(x) x[[1]]))
sour_pos = grep("TDC", sour)
if(sum(sour[sour_pos] != pas.list) == 0){
  sour[sour_pos] = pas
}else{
  print("There are something wrong!")
}
write.table(sour, file = "./reviced_NCAM_flipped_MAF0.05_chr22_rehh.inp", quote = F, col.names = F, row.names = F)
