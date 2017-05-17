## IN ALPS
setwd("~/TaiwanBiobank/NGS_illumina_iSAAC_vcf")
# Getting individual list
# lftp cloud-ftp -e "set ftp:passive-mode true;cd /TaiwanBioBank/NGS_illumina_iSAAC_vcf/iSAAC_VCF;ls>list.txt;bye"
new_list = read.table("./list.txt", stringsAsFactors = FALSE)
temp = new_list$V9
tem_tb = data.frame(temp)
write.table(tem_tb, file = "./ind_list_for_merge.txt", quote = F, col.names = F, row.names = F)

## Turn vcf to PLINK bfiles
# record_sheet = data.frame()
for(i in 1:100){
    lftp = paste0("lftp cloud-ftp -e \"set ftp:passive-mode true;get ~/TaiwanBioBank/NGS_illumina_iSAAC_vcf/iSAAC_VCF/", temp[i],"/",temp[i],"_S1.vcf",";bye\"")
    system(lftp)
    system(paste0("plink -vcf ", temp[i],"_S1.vcf --snps-only --exclude dot_rm.txt --double-id --make-bed --out ", temp[i]))
    system(paste0("rm ", temp[i],"_S1.vcf"))
    system(paste0("mkdir ",temp[i],"_bfile"))
    system(paste0("mv ",temp[i],".* /home/u00ccy01/TaiwanBiobank/NGS_illumina_iSAAC_vcf/",temp[i],"_bfile"))
    lftp = paste0("lftp cloud-ftp -e \"set ftp:passive-mode true;set ssl:verify-certificate no;mirror /home/u00ccy01/TaiwanBiobank/NGS_illumina_iSAAC_vcf/",temp[i],"_bfile -R /TaiwanBioBank/NGS_illumina_iSAAC_vcf/bfile;bye\"")    
    system(lftp)
    system(paste0("rm -R ",temp[i],"_bfile"))
    # start_time = Sys.time()
    # end_time = Sys.time()
    # duration = end_time - start_time
    # record_sheet[i,1] = ind_list[i,]
    # record_sheet[i,2] = format(start_time, "%H:%M:%S")
    # record_sheet[i,3] = format(end_time, "%H:%M:%S")
    # record_sheet[i,4] =  round(duration)
    #file.remove()
    # if(nrow(record_sheet) == 2){
    #     write.table(record_sheet, file = "./download_time.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
    # }
}

## Check error
## IN ALPS
# lftp cloud-ftp -e "set ftp:passive-mode true;cd /TaiwanBioBank/NGS_illumina_iSAAC_vcf/bfile;ls>list.txt;bye"
## IN Mac
setwd("~/TaiwanBioBank/")
bfile_list = read.table("./bfile_list.txt", stringsAsFactors = FALSE)
head(bfile_list)
fam = bfile_list$V9[grep(".fam", bfile_list$V9)]
length(fam)
head(fam)
library(stringr)
fam = str_sub(fam, 1, 12)
log = bfile_list$V9[grep(".log", bfile_list$V9)]
length(log)
head(log)
log = str_sub(log, 1, 12)
error = setdiff(log, fam)
# Errors of "NGS20140705F" and "NGS20140711G" were caused by missing vcf file
# Re-converting
which(new_list$V9 %in% "NGS20140705F")
which(new_list$V9 %in% "NGS20140711G")
for(i in c(123, 168)){
    lftp = paste0("lftp cloud-ftp -e \"set ftp:passive-mode true;get ~/TaiwanBioBank/NGS_illumina_iSAAC_vcf/iSAAC_VCF/", temp[i],"/",temp[i],"_S1.vcf",";bye\"")
    system(lftp)
    system(paste0("plink -vcf ", temp[i],"_S1.vcf --snps-only --exclude dot_rm.txt --double-id --make-bed --out ", temp[i]))
    system(paste0("rm ", temp[i],"_S1.vcf"))
    system(paste0("mkdir ",temp[i],"_bfile"))
    system(paste0("mv ",temp[i],".* /home/u00ccy01/TaiwanBiobank/NGS_illumina_iSAAC_vcf/",temp[i],"_bfile"))
    lftp = paste0("lftp cloud-ftp -e \"set ftp:passive-mode true;set ssl:verify-certificate no;mirror /home/u00ccy01/TaiwanBiobank/NGS_illumina_iSAAC_vcf/",temp[i],"_bfile -R /TaiwanBioBank/NGS_illumina_iSAAC_vcf/bfile;bye\"")    
    system(lftp)
    system(paste0("rm -R ",temp[i],"_bfile"))
}
# Re-writing individuals list
# IN MAC
error = error[!(error %in% c("NGS20140705F", "NGS20140711G"))]
error = sub("\\.+", "", error)
ind_list = read.table("./ind_list_for_merge.txt", stringsAsFactors = F)
ind_list = ind_list$V1[-(which(ind_list$V1 %in% error))]
ind_list_df = data.frame(ind_list)
write.table(ind_list_df, file = "./new_merge_list.txt", quote = F, row.names = F, col.names = F)


# install.packages("RCurl")
# library(RCurl)
# un = readline("Type the username:")
# pw = readline("Type the password:")
# upw = paste(un, pw, sep = ":")
# getURL(url_list[1], userpwd = upw)
# 
# download.file("ftp://alps.st.nchc.tw:2122/TaiwanBioBank/NGS_illumina_iSAAC_vcf/iSAAC_VCF/NGS20140601B/NGS20140601B_S1.genome.vcf.gz", destfile = "./test.vcf.gz")
# 
# getURL(URLencode("ftp://alps.st.nchc.tw:2122/TaiwanBioBank/NGS_illumina_iSAAC_vcf/iSAAC_VCF/NGS2015032H/NGS2015032H_S1.genome.vcf.gz"), userpwd = "u00ccy01:yo@4100030307")
# ftpUpload("./list.txt", "ftp://chunyo:yo80106@120.126.38.111/Public/chunyo_test/list.txt")