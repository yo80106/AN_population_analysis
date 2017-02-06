### Per-individual QC
## Data preprocessing
system("plink --file 96samples --make-bed --out 96samples")
system("plink --file 8samples --make-bed --out 8samples")
system("plink --bfile 96samples --remove list.txt --make-bed --out extract")
system("plink --bfile 8samples --bmerge extract.bed extract.bim extract.fam --recode --out taitung")

## Sex check
fam = read.table("./new.fam")
sex = read.csv("./96samples_sex_list.csv",header = F)
# To find out the ID names whether consist between two files or not 
fam$V5 = sex$V3
fam$V2 == sex$V2
# Write out new fam file which record the right sex information
levels(fam$V5)[1] = 2
levels(fam$V5)[2] = 1
write.table(fam,"./new.fam",col.names = F,row.names = F,quote = F)
sex_check("taitung", "taitung")

## Missing data rates and outlying heterozygosity
missing_het_ind("taitung", "96samples_sex_list.csv", "taitung")

## IBD
IBD("taitung","taitung")

## Removal of all individuals failing QC
ind_qc_rm("taitung", "taitung")

## To remove all the individuals failed in the previous QC steps.
system("plink --bfile taitung --remove taitung_fail_ind_qc.txt --make-bed --out clean_inds_data")

### Per-marker QC
## Excessive missing data rate
missing_snp("clean_inds_data", "clean_inds_data")

## Plot HWE p-value distribution
hwe_test("clean_inds_data", "clean_inds_data")

# Pie chart
# pdf("pie_missing_geno.pdf")
lmiss = read.table("clean_inds_data.lmiss", header = T)
small_table = data.frame(c(sum(lmiss$F_MISS == 0),sum(lmiss$F_MISS < 0.03)-
            sum(lmiss$F_MISS == 0),sum(lmiss$F_MISS>0.03))
            ,row.names = c("missingness = 0","0 < missingness < 0.03","missingness > 0.03"))
colnames(small_table) = "counts"


lbls = paste(row.names(small_table),"\n",small_table$counts,sep = "")
pie(small_table$counts,labels = lbls,col = rainbow(length(lbls)))
# dev.off()

## Removal of all markers failing QC
## PLINK SECTION
system("plink --bfile clean_inds_data --maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out clean_inds_snp_data")