setwd("~")
library(stringr)
pan_asia_fam = read.table("~/PCA/R_analysis/HGDP_Pan_merge/Genotypes_All.fam")
pan_asia_fam$pop = NA

for(i in 1:nrow(pan_asia_fam)){
    pop =  str_sub(pan_asia_fam$V2[i],start = 1,end = 5)
    pan_asia_fam$pop[i] = pop 
}
which(pan_asia_fam$pop == "CEU-N")
pan_asia_fam$pop[1720:nrow(pan_asia_fam)] = 
    str_sub(pan_asia_fam$pop[1720:nrow(pan_asia_fam)],start = 1, end = 3)
unique(pan_asia_fam$pop)

setwd("~/PCA/R_analysis/location/")
for(i in unique(pan_asia_fam$pop)){
    url = paste("http://www4a.biotec.or.th/export_pasnp/Pop_Genotype/",i,"_information.txt",sep="")
    file_name = paste(i,".txt",sep="")
    print(file_name)
    download.file(url,file_name)
}
file_list = list.files("./",pattern = ".txt")
for(file in file_list){
    if(!exists("location_list")){
        location_list = read.delim(file,stringsAsFactors = FALSE)
    }
    if(exists("location_list")){
        temp_list = read.delim(file,stringsAsFactors = FALSE)
        location_list = rbind(location_list,temp_list)
        rm(temp_list)
    }
}
lat.long = location_list[,c(1,7,8)]
lat.long$lat.decimal = NA
lat.long$long.decimal = NA
convert = function(vector){
    (vector[1] + (vector[2]/60) + (vector[3]/3600))*vector[4]
}
count = 1
for(lat in lat.long$Latitude){
    lat = gsub("deg","",lat)
    lat = gsub("'","",lat)
    lat = gsub("N",1,lat)
    lat = gsub("S",-1,lat)
    lat.num = as.numeric(strsplit(lat," ")[[1]])
    lat.long$lat.decimal[count] = convert(lat.num)
    count = count + 1
}
count = 1
for(lon in lat.long$Longitude){
    lon = gsub("deg","",lon)
    lon = gsub("'","",lon)
    lon = gsub("E",1,lon)
    lon = gsub("W",-1,lon)
    lon.num = as.numeric(strsplit(lon," ")[[1]])
    lat.long$long.decimal[count] = convert(lon.num)
    count = count + 1
}
lat.long = lat.long[!(lat.long$Pop_name %in% c("CEU","JPT","YRI")),]

hgdp_location = read.table("~/PCA/R_analysis/HGDP_localist.txt",stringsAsFactors = FALSE, header = TRUE)
colnames(hgdp_location) = c("Population","Region","Latitude","Longitude")
hgdp_location = hgdp_location[!(hgdp_location$Population %in% c("Han","Han−NChina","Miao","Japanese")),]
# names(lat.long)[c(1,4,5)] = names(hgdp_location)
# lat.long = rbind(hgdp_location,lat.long[,c(1,4,5)])

## Making geographic distribution of sampling sites
setwd("~/PCA/R_analysis/")
install.packages("maps")
install.packages("mapdata")
library(maps)
library(mapdata)

pdf(file = "~/PCA/R_analysis/map.pdf")
map('worldHires')
points(lat.long$long.decimal,lat.long$lat.decimal,col = "blue",pch = 20,cex = 0.6)
#text(lat.long$long.decimal+4,lat.long$lat.decimal, labels = lat.long$Pop_name,cex = 0.2,col = "blue")
points(hgdp_location$Longitude,hgdp_location$Latitude,col = "red",pch = 20,cex = 0.6)
#text(hgdp_location$Longitude+3,hgdp_location$Latitude, labels = hgdp_location$Population,cex = 0.2,col = "red")
dev.off()