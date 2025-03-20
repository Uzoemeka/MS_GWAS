#________________________________________________
#____________ Extract Top SNP from GWAS
#________________________________________________

# Load GWAS results for 10-year and 15-year predictions
gwas_10 <- fread("/gwas_results/No_year_and_dmt_pheno.pred_10yr_edss.assoc.linear")
gwas_15 <- fread("/Karim/gwas_results/No_year_and_dmt_pheno.pred_15yr_edss.assoc.linear")

# Sort by p-value and extract top 200 SNPs
gwas_10 <- gwas_10[order(gwas_10$P),]
gwas_15 <- gwas_15[order(gwas_15$P),]

# Select top 200 SNPs
gwas_10_1 <- gwas_10[1:200,]
gwas_15_1 <- gwas_15[1:200,]

# Write top SNPs to text files
write.table(gwas_10_1$SNP, "/Karim/gwas_results/No_year_and_dmt_Pred_10yr_top_200.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(gwas_15_1$SNP, "/Karim/gwas_results/No_year_and_dmt_Pred_15yr_top_200.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Compare top SNPs between 10-year and 15-year predictions
length(intersect(gwas_15_1$SNP, gwas_10_1$SNP))
length(setdiff(gwas_15_1$SNP, gwas_10_1$SNP))

# Combine and remove duplicates
pred_10_15 <- rbind(gwas_15_1[,1:2], gwas_10_1[,1:2])
length(unique(pred_10_15$SNP))

# Write combined top SNPs to a new file
write.table(unique(pred_10_15$SNP), "/MS_GWAS/Karim/gwas_results/No_year_and_dmt_Pred_10yr_15yr_top.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#________________________________________________
# Extract Genotype Data for Selected SNPs
#________________________________________________

# Define input files
input_file <- "/MS_GWAS/Karim/ped_map_files/karim_clean_threemonths_imputed.raw"
columns_file <- "/MS_GWAS/Karim/gwas_results/No_year_and_dmt_Pred_10yr_15yr_top.txt"

# Step 1: Extract the header line and convert to a newline-separated list
system(paste("head -n 1", input_file, "| tr ' ' '\n' > /MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_header_columns.txt"))

# Step 2: Process the SNPs for genotype extraction
pheno_column <- fread(columns_file, header = FALSE)
geno_column <- fread("/MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_header_columns.txt", header = FALSE)

# Extract matching SNPs from the header
geno_column$new <- sub("_.*", "", geno_column$V1)
geno_column_new <- geno_column[geno_column$new %in% pheno_column$V1,]
write.table(geno_column_new$V1, "/MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_new.txt", quote = F, col.names = F, row.names = F, sep = "\t")

# Step 3: Find the indices of the desired columns
columns_file_new <- "/MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_new.txt"
system(paste("grep -n -f", columns_file_new, "/MLM/MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_header_columns.txt", "| cut -d: -f1 > /Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_column_indices.txt"))

# Step 4: Extract the desired columns
indices <- paste(scan("MS_GWAS/Karim/No_year_and_dmt_Pred_10yr_15yr_top_column_indices.txt", what = "", sep = "\n"), collapse = ",")
system(paste("cut -d ' ' -f", indices, input_file, ">", "/MS_GWAS/Karim/ped_map_files/No_year_and_dmt_top_10_15yr_imputed.raw"))

# Optional: Extract specific columns using awk for large files
# Example: Extracting the first 50,000 columns
system(paste("awk '{for(i=1; i<=50000 && i<=NF; i++) printf \"%s \", $i; print \"\"}'", input_file, ">", "/MLM_GWAS/karim_clean_threemonths_imputed_noFID_split_1.raw"))

# Example: Extracting the first column and columns in between (for large files)
system(paste("awk '{printf \"%s\", $1; for (i = 5100001; i <= 5117076; i++) printf \"\t%s\", $i; printf \"\\n\"}'", input_file, ">", "/MLM_GWAS/karim_clean_threemonths_imputed_noFID_split_103.raw"))
