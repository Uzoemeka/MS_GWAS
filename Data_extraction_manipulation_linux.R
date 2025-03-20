#_______________________________________________________________________________
#     Converting and extracting relevant data
#_______________________________________________________________________________

# Extract SNP from bed, bim, fam file
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--bfile /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/karim_clean \
--extract /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/snp_ids.txt \
--make-bed --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/karim_rs10191329")

# Convert bed, bim, fam file to ped and map files
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--bfile /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/karim_rs10191329 \
--recode --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_rs10191329")

# Extract SNP and convert to map and ped files
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--file /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean \
--extract /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/snp_ids.txt \
--recode --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_rs10191329")

# Convert ped and map file to .raw file containing 0,1,2
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--file /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_rs10191329 \
--recodeA --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_rs10191329_3")

# Convert SNP from bed, bim, fam file to map and ped file
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--bfile /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/karim_clean_threemonths \
--recode --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths")

# Convert map and ped file to vcf file
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--file /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths \
--recode vcf --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths")

# Run Beagle to impute missing values
system("java -jar /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/beagle.21Apr21.304.jar \
gt=/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths.vcf \
out=/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths_imputed")

# Convert imputed VCF back to ped/map format using plink
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--vcf /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths_imputed.vcf.gz \
--recode --const-fid 0 \
--out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths_imputed")

# Convert SNP from map and ped file to .raw file containing 0,1,2
system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--file /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths \
--recodeA --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths")

system("/Users/uzoemeka/Documents/plink_mac_20231211/plink \
--file /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths_imputed \
--recodeA --out /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths_imputed")
