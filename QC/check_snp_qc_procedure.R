setwd("~")
# Equal to --bfile new --maf 0.01 --make-bed --out maf_control
no_clean_frq = read.table("~/QC/test/snp_qc_all/maf_test.frq",header = T)
table(is.na(no_clean_frq$MAF))
no_total_fail_MAF = sum(no_clean_frq$MAF < 0.01, na.rm = TRUE)
# Equal to --bfile new --hwe 0.00001 --make-bed --out hew_control
no_clean_hwe = read.table("~/QC/test/snp_qc_all/hwe_test.hwe",header = T)
table(is.na(no_clean_hwe$P))
no_total_fail_hwe = sum(no_clean_hwe$P<0.00001, na.rm = TRUE)
# Equal to --bfile new --geno 0.03 --make-bed --out geno_control
no_clean_geno = read.table("~/QC/test/snp_qc_all/geno_test.lmiss",header = T)
table(is.na(no_clean_geno$F_MISS))
no_total_fail_geno = sum(no_clean_geno$F_MISS>=0.03)

### Equal to --bfile new --geno 0.03 --hwe 0.00001 
### --maf 0.01 --make-bed --out control
## In this case, 
## 70144 variants removed due to missing genotype data (--geno)
## 1277888 variants removed due to MAF threshold(s) (--maf)
## 112 variants removed due to Hardy-Weinberg exact test (--hwe)

# PLINK will first execute --geno command to remove high missing genotype
# Results are same as above.
sum(no_clean_geno$F_MISS>=0.03)
# PLINK then execute --maf after doing --geno.
geno_control = read.table("~/QC/test/snp_qc_all/test/maf_control.frq"
                          , header = TRUE)
table(is.na(geno_control$MAF))
no_total_fail_MAF = sum(geno_control$MAF < 0.01, na.rm = TRUE)
no_total_fail_MAF
# PLINK finally execute --hwe after doing --geno
hwe_control = read.table("~/QC/test/snp_qc_all/test/hwe_control.hwe"
                         , header = TRUE)
table(is.na(hwe_control$P))
no_total_fail_hwe = sum(hwe_control$P < 0.00001, na.rm = TRUE)