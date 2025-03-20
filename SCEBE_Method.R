#_______________________________________________________________________________
#                             SCEBE - Two-stage Multilevel Model GWAS
#_______________________________________________________________________________

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(MASS)
library(nlme)
library(Matrix)
library(lme4)
library(memoise)
library(dplyr)  # To revive the isFALSE() function for sim_slopes()
library(tidyr)
library(stringr)
library(SCEBE)  # Custom package for SCEBE method

# Prepare the data
colnames(papisQ_30)[1] <- "ID"
final_data <- papisQ_30[complete.cases(papisQ_30),]  # Remove rows with missing values

# Check the structure of the cleaned data
length(unique(final_data$ID))
head(final_data)

# Fit the mixed-effects model (LME)
m <- lmer(edssScore ~ 1 + time + (1 + time | ID), 
          data = final_data, REML = FALSE, 
          control = lmerControl(optimizer = "bobyqa"))
print(m)

# Prepare phenotype and genotype data for SCEBE analysis
dat <- final_data[, c("ID", "edssScore", "time", "s.time", "l.time", "sex.cat", "DMT.cat", 
                      "ageonset.cent", "PC1", "PC2", "PC3", "PC4", "PC5")]
head(dat)

genodat <- unique(final_data[, c(1, 14:364)], by = "ID")  # Genotype data
head(genodat[, 1:5], 2)

# Run SCEBE analysis
sm <- scebe_sim(phenoData = dat, genoData = genodat[,-1], m, 
                Time = "time", pheno = "edssScore", method = "scebe")
head(sm)

#_______________________________________________________________________________
#                            SCEBE Results Processing
#_______________________________________________________________________________

# Load SCEBE GWAS results
paths_twas <- fs::dir_ls("/Karim/gwas_results")
paths_twas <- paths_twas[grepl("scebe_random_slope_gwas.txt", paths_twas)]

# Read the TWAS results
data_twas <- map(paths_twas, read.table, header = TRUE)
data_twas <- bind_rows(data_twas)  # Combine all data into a single dataframe
dim(data_twas)

# Process TWAS results
data_twas$chr <- as.numeric(gsub("chr", "", sub("_.*$", "", rownames(data_twas))))
data_twas$pos <- as.numeric(sub("_.*$", "", sub("^([^_]*_){1}", "", rownames(data_twas))))
data_twas$A1 <- sub("_.*$", "", sub("^([^_]*_){2}", "", rownames(data_twas)))
data_twas$A2 <- sub("_.*$", "", sub("^([^_]*_){3}", "", rownames(data_twas)))
data_twas$SNP <- paste(data_twas$chr, data_twas$pos, sep = ":")
data_twas$N <- rep(905, dim(data_twas)[1])  # Sample size
head(data_twas)

# Load and process the BIM file for SNP information
library(data.table)
bim_file <- fread("/karim_clean_threemonths_snpids.updated38.bim")
bim_file$SNP <- paste(bim_file$V1, bim_file$V4, sep = ":")
head(bim_file)

# Merge TWAS results with SNP data
data_twas <- merge(data_twas, bim_file, by = "SNP")
head(data_twas)

# Extract and order data for intercept (BETA1)
data_intercept <- data_twas[, c('chr', 'pos', 'V2', 'A1', 'A2', 'est.sc.b1', 'se.sc.b1', 'p.sc.b1', 'N')]
colnames(data_intercept) <- c('CHR', 'POS', "SNP", 'A1', 'A2', 'BETA', 'SE', 'P', 'N')
data_intercept <- data_intercept[order(data_intercept$P),]
head(data_intercept)

# Extract and order data for slope (BETA2)
data_slope <- data_twas[, c('chr', 'pos', 'V2', 'A1', 'A2', 'est.sc.b2', 'se.sc.b2', 'p.sc.b2', 'N')]
colnames(data_slope) <- c('CHR', 'POS', "SNP", 'A1', 'A2', 'BETA', 'SE', 'P', 'N')
data_slope <- data_slope[order(data_slope$P),]
head(data_slope)

#_______________________________________________________________________________
#                            Plotting Results
#_______________________________________________________________________________

# Load QQMAN package for Manhattan plot
library(qqman)

# Plot Manhattan plot for random slope
ggsave(plot = manhattan(data_slope, chr = "CHR", bp = "POS", p = "P", snp = "SNP", 
                        col = c("skyblue", "tomato"), suggestiveline = -log10(1e-5), 
                        genomewideline = -log10(5e-8), main = "SCEBE Random slope", 
                        annotateTop = TRUE, ylim = c(0, 10), cex = 0.6), 
       "manhattan_plotb.pdf", width = 10, height = 6)

# QQ plot for p-values of the random slope
qq(data_slope$P, main = "SCEBE Random slope")

#_______________________________________________________________________________
#                            Comparing GWAS Results
#_______________________________________________________________________________

# Load MLM GWAS results
mlm_gwas <- fread("/gwas_results/MLM_direct_gwas_results_top_hits_imputed.txt")
head(mlm_gwas)

# Load PLINK GWAS results (10-year and 15-year)
plink_gwas10 <- fread("/gwas_results/No_year_and_dmt_pheno.pred_10yr_edss.assoc.linear")
plink_gwas15 <- fread("/gwas_results/No_year_and_dmt_pheno.pred_15yr_edss.assoc.linear")

# Process PLINK GWAS results (10-year and 15-year)
plink_gwas10 <- plink_gwas10[order(plink_gwas10$P),]
plink_gwas15 <- plink_gwas15[order(plink_gwas15$P),]

# Select top 200 SNPs
plink_gwas10_1 <- plink_gwas10[1:200,]
plink_gwas15_1 <- plink_gwas15[1:200,]

# Create unique variable identifiers for merging
plink_gwas10_1$variable <- paste(plink_gwas10_1$SNP, ":", plink_gwas10_1$A1, sep = "")
plink_gwas15_1$variable <- paste(plink_gwas15_1$SNP, ":", plink_gwas15_1$A1, sep = "")

# Remove ":" and replace with "_"
plink_gwas10_1$variable <- gsub(":", "_", plink_gwas10_1$variable)
plink_gwas15_1$variable <- gsub(":", "_", plink_gwas15_1$variable)

# Merge MLM with PLINK results for 10 and 15 years
mlm_gwas_plink_gwas10 <- merge(mlm_gwas, plink_gwas10_1, by = "variable")
mlm_gwas_plink_gwas15 <- merge(mlm_gwas, plink_gwas15_1, by = "variable")

# Plot P-value comparisons (LME vs Standard GWAS)
ggplot(mlm_gwas_plink_gwas10, aes(x = -log10(p_value_10), y = -log10(P))) +
  geom_point(size = 3, alpha = 0.6) +
  labs(title = "Scatter plot of P-values from LME vs standard GWAS (10 years post-onset)",
       x = "-log10(Pvalue) LME", y = "-log10(Pvalue) standard GWAS") +
  theme_minimal()

# Plot effect size comparisons (LME vs Standard GWAS)
ggplot(mlm_gwas_plink_gwas10, aes(x = coefficient_10, y = BETA)) +
  geom_point(size = 3, alpha = 0.6) +
  labs(title = "Scatter plot of Effect sizes from LME vs standard GWAS (10 years post-onset)",
       x = "Effect sizes LME", y = "Effect sizes standard GWAS") +
  theme_minimal()

# Plot P-value comparisons (15 years post-onset)
ggplot(mlm_gwas_plink_gwas15, aes(x = -log10(p_value_15), y = -log10(P))) +
  geom_point(size = 3, alpha = 0.6) +
  labs(title = "Scatter plot of P-values from LME vs standard GWAS (15 years post-onset)",
       x = "-log10(Pvalue) LME", y = "-log10(Pvalue) standard GWAS") +
  theme_minimal()

