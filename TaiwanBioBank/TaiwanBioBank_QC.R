setwd("~/TaiwanBioBank/raw_data/bfiles/Archive/")
source("~/Git/QC/QC_functions.R")

## Ind QC
sex_check("TWBR10505-02_update", "TWBB")

TW_fam = read.table("./TWBR10505-02_update.fam", stringsAsFactors = FALSE)
head(TW_fam)
TW_fam$pop = "Han"
write.table(TW_fam[,c(1,2,7)], file = "./TW_pop.txt", quote = F, row.names = F, col.names = F)
missing_het_ind("TWBR10505-02_update", "TW_pop", "TWBB")

IBD("TWBR10505-02_update", "TWBB")

ind_qc_rm("TWBB", "TWBB")
system("plink --bfile TWBR10505-02_update --remove TWBB_fail_ind_qc.txt --make-bed --out TWBB_rm_ind")

ind_qc_info("TWBB", "TWBB")

## SNP QC
missing_snp("TWBR10505-02_update", "TWBB")
TW_lmiss = read.table("./TWBB.lmiss", stringsAsFactors = FALSE, header = TRUE)
head(TW_lmiss)
sum(TW_lmiss$F_MISS > 0.03)

hwe_test("TWBR10505-02_update", "TWBB")
TW_hwe = read.table("./TWBB.hwe", stringsAsFactors = FALSE, header = TRUE)
system("plink --bfile TWBB_rm_ind --maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out TWBB_rm_ind_snp")
