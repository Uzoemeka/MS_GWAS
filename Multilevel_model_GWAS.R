# =============================================================================
#                Multilevel Model for Genotypes
# =============================================================================

# Loading necessary libraries
library(tidyr)
library(lattice)
library(lmtest)
library(readxl)
library(car)
library(backports)     # For reviving the isFALSE() function for sim_slopes()
library(effects)
library(ggplot2)
library(interactions)
library(lme4)
library(lmerTest)
library(psych)
library(plyr)
library(dplyr)
library(R2MLwiN)
library(data.table)
library(writexl)
library(haven)
library(tidyverse)
library(nlme)
library(curl)
library(rstan)
library(brms)
library(lubridate)
library(kableExtra)
library(foreign)
library(multcomp)
library(groupdata2)
library(merTools)

# Setting MLwiN path
options(MLwiN_path = '/path/mlnscript')

# =============================================================================
#                             Data Preparation
# =============================================================================

# Load dataset
papis <- read_dta("Edssdatatata2nd_three_month.dta")
papis <- as.data.frame(papis)

# Renaming column
names(papis)[names(papis) == 'date'] <- 'relapseDate'

# Subsetting and reordering columns
papis <- papis[, c("patientIdentifier", "dateTime", "msDateDiagnosis", "msDateOnset", 
                   "DMT", "relapseDate", "edssScore", "edssScoreDecimal", "msEventType", 
                   "sex", "dob")]

# Renaming columns for clarity
papis$patid <- papis$patientIdentifier

# Converting date columns to Date format
papis$msDateDiagnosis <- as.Date(papis$msDateDiagnosis, "%d/%m/%Y")
papis$msDateOnset <- as.Date(papis$msDateOnset, "%d/%m/%Y")
papis$dateTime <- as.Date(papis$dateTime, "%d/%m/%Y")
papis$dob <- as.Date(papis$dob, "%d/%m/%Y")
papis$relapseDate <- as.Date(papis$relapseDate, "%d/%m/%Y")

# =============================================================================
#                      Time and Age Calculations
# =============================================================================

# Create time-related variables
papis$ageonset <- as.numeric((papis$msDateOnset - papis$dob) / 365.25)
papis$timeage <- as.numeric((papis$dateTime - papis$dob) / 365.25)
papis$timediag <- as.numeric((papis$dateTime - papis$msDateDiagnosis) / 365.25)
papis$timeonset <- as.numeric((papis$dateTime - papis$msDateOnset) / 365.25)

# Filter out negative timeonset values
papis <- papis[papis$timeonset >= 0, ]

# Creating diagnosis year and categorical year variables
papis <- papis %>% mutate(year_diagnosis = as.numeric(format(msDateDiagnosis, '%Y')))
papis$year_diagnosis_cat <- as.factor(
  ifelse(papis$msDateDiagnosis > as.POSIXct("1953", format = "%Y") & papis$msDateDiagnosis < as.POSIXct("2005", format = "%Y"), "year1",
         ifelse(papis$msDateDiagnosis > as.POSIXct("2006", format = "%Y") & papis$msDateDiagnosis < as.POSIXct("2010", format = "%Y"), "year2", "year3"))
)

# =============================================================================
#                      Creating Nonlinear Time Terms
# =============================================================================

# Nonlinear time transformations
papis$time <- papis$timeonset + 1
papis$l.time <- log(papis$time)
papis$s.time <- sqrt(papis$time)

# Factorizing and releveling categorical variables
papis$DMT <- as.factor(papis$DMT)
papis$DMT <- relevel(papis$DMT, ref = "notDMT")

papis$sex <- as.factor(papis$sex)
papis$sex <- relevel(papis$sex, ref = "MALE")

# Create categorical variables for sex and DMT
papis$sex.cat <- ifelse(papis$sex == "MALE", 0, 1)
papis$DMT.cat <- ifelse(papis$DMT == "notDMT", 0, 1)

# Create categorical variable for year of diagnosis
papis$year.diagnosis.cat <- ifelse(papis$year_diagnosis_cat == "year1", 0,
                                   ifelse(papis$year_diagnosis_cat == "year2", 1, 2))

# =============================================================================
#                 Data Aggregation by Quarter
# =============================================================================

# Aggregating observations within a quarter-year period
id <- unique(papis$patid)
papis_quarter <- NULL

for (i in id) {
  papis_i <- as.data.frame(papis %>% filter(patid == i) %>%
                             select(patid, dateTime, edssScore, ageonset, timeage, timeonset, time, s.time, l.time, sex.cat, DMT.cat, year.diagnosis.cat))
  
  break_dates <- seq.Date(min(papis$dateTime), as.Date("2022-01-01"), by = "quarter")
  cut_i <- cut(papis_i$dateTime, breaks = break_dates)
  papis_median <- aggregate(papis_i, by = list(cut_i), median)
  
  papis_quarter <- rbind(papis_quarter, papis_median)
}

# Correcting DMT.cat values
papis_quarter$DMT.cat[papis_quarter$DMT.cat == 0.5] <- 1

# Check number of unique patients
length(unique(papis_quarter$patid))

# =============================================================================
#                         PCA Data
# =============================================================================

# Load PCA data and format
pca <- fread("karim_clean_threemonths_pca1.eigenvec")
pca <- pca[, 1:7]
pca$ids <- sub("^[^_]*_", "", pca$FID)

# Check first two rows and number of unique IDs
head(pca, 2)
length(unique(pca$FID))



# =============================================================================
#                   Data Preparation and Merging
# =============================================================================

# Set directory path for the raw data
di <- "Raw_data_split_3"

# Load raw genotype data
raw_file <- fread(paste(di, "/karim_clean_threemonths_imputed_noFID_split_103.raw", sep = ""))

# Display first 2 rows and dimensions of raw_file
head(raw_file[, 1:10], 2)
dim(raw_file)

# Extract patient IDs from IID column
raw_file$ids <- sub("^[^_]*_", "", sub("^[^_]*_", "", sub("^[^_]*_", "", sub("^[^_]*_", "", raw_file$IID))))

# Check unique IDs in raw_file
length(unique(raw_file$IID))
length(unique(raw_file$ids))

# Load all genotyped patient IDs
ids <- fread("all.ids_Genotyped_patients_threemonthspostrelapse.csv")
ids$ids <- ids$IID

# Merge raw genotype data with patient IDs
rawfile <- merge(raw_file, ids, by = "ids")

# Check unique IDs after merge
length(unique(rawfile$ids))

# Remove unwanted columns (duplicate or unnecessary identifiers)
rawfile <- rawfile[,-c("IID.x", "IID.y", "sex.cat", "ageonset", "ageonset.cent", "dupl.x", "FID", "PatID", "MatID", "SEX", "PHENO", "dupl.y")]

# Standardize the data (scale all columns except "patid" and "ids")
rawfile <- rawfile %>% mutate(across(.cols = -all_of(c("patid", "ids")), .fns = scale))

# Rename columns by replacing ':' with '_'
colnames(rawfile)[-c(1, dim(rawfile)[2])] <- gsub(":", "_", colnames(rawfile)[-c(1, dim(rawfile)[2])])

# Display first 3 rows of the scaled data
head(rawfile[, 1:10], 3)
dim(rawfile)[2]

# Check number of unique patient IDs in rawfile
length(unique(rawfile$ids))

# =============================================================================
#                   Merge PCA Data and Further Processing
# =============================================================================

# Merge rawfile with PCA data based on 'ids'
rawfile1 <- merge(rawfile, pca, by = "ids")

# Remove PCA-related columns after merging
rawfile1 <- rawfile1[,-c("ids", "IID", "FID")]

# Display first 2 rows and dimensions of rawfile1
head(rawfile1[, 1:10], 2)
dim(rawfile1)

# Check number of unique patient IDs in rawfile1
length(unique(rawfile1$patid))

# =============================================================================
#                   Merge with Papis_quarter Data
# =============================================================================

# Merge rawfile1 with papis_quarter data based on 'patid'
papis_quarter <- merge(papis_quarter, rawfile1, by = "patid")

# Display first 2 rows and dimensions of papis_quarter
head(papis_quarter[, 1:10], 2)
dim(papis_quarter)

# Check number of unique patient IDs in papis_quarter
length(unique(papis_quarter$patid))


# =============================================================================
#                      Defining Other Covariates
# =============================================================================

# Convert categorical variables to factors with labels
papis_quarter$sex.cat <- factor(papis_quarter$sex.cat, labels = c("male", "female"))
papis_quarter$DMT.cat <- factor(papis_quarter$DMT.cat, labels = c("offDMT", "onDMT"))
papis_quarter$year.diagnosis.cat <- factor(papis_quarter$year.diagnosis.cat, labels = c("year1", "year2", "year3"))
papis_quarter$yeardiagnosiscat <- factor(papis_quarter$year.diagnosis.cat, labels = c("year1", "year2", "year3"))

# Centering ageonset by subtracting the mean (31.5)
papis_quarter$ageonset.cent <- papis_quarter$ageonset - 31.5

# Final dataset (no missing data)
final_data <- papis_quarter  # Can uncomment to filter rows with complete cases: [complete.cases(papis_quarter),]
length(unique(final_data$patid))  # Check number of unique patient IDs

# =============================================================================
#                    30-Year Truncated Data
# =============================================================================

# Filter data for patients with time <= 30 years
papisQ_30 <- as.data.frame(final_data %>% group_by(patid) %>% filter(time <= 30))
length(unique(papisQ_30$patid))  # Check number of unique patient IDs in 30-year data

# Create new variables by copying from existing ones for easier analysis
papisQ_30$sexcat <- papisQ_30$sex.cat
papisQ_30$DMTcat <- papisQ_30$DMT.cat
papisQ_30$stime <- papisQ_30$s.time
papisQ_30$ltime <- papisQ_30$l.time
papisQ_30$ageonsetcent <- papisQ_30$ageonset.cent
papisQ_30$yeardiagnosiscat <- papisQ_30$year.diagnosis.cat

# =============================================================================
#                    Selecting Variables for Analysis
# =============================================================================

# Exclude unnecessary variables for modeling or analysis
variables <- names(papisQ_30)[!names(papisQ_30) %in% c(
  "Group.1", "ageonset", "timeage", "timeonset", "time", "s.time", "l.time", 
  "dateTime", "sex.cat", "patid", "edssScore", "ageonset.cent", "year.diagnosis.cat", 
  "DMT.cat", "yeardiagnosiscat", "sexcat", "DMTcat", "stime", "ltime", 
  "ageonsetcent", "PC1", "PC2", "PC3", "PC4", "PC5"
)]

# =============================================================================
#                      Load Necessary Libraries
# =============================================================================

# Parallel computing, mixed models, and tidying model output
library(parallel)  # For parallel computing
library(lme4)  # For linear mixed-effects models
library(broom.mixed)  # For tidying mixed-model outputs

# =============================================================================
#              Define Variables and Model Setup
# =============================================================================

# Constant variables for the model formula
constant_vars <- c(
  "1", "s.time*ageonset.cent", "l.time*ageonset.cent", 
  "s.time*sex.cat", "l.time*sex.cat", "PC1", "PC2", "PC3", "PC4", "PC5"
)

# Variables to loop over in the model
loop_vars <- variables  # List of looped variables

# Random effects specification for the model
random_effect <- "(1 + s.time + l.time | patid)"  # Random intercept and slopes for `s.time` and `l.time`

# =============================================================================
#              Define Function to Fit the Model and Extract Results
# =============================================================================

fit_model <- function(var) {
  # Create the model formula dynamically by adding interaction terms with the loop variable
  formula <- as.formula(paste(
    "edssScore ~", 
    paste(c(constant_vars, paste("l.time", "*", var), paste("s.time", "*", var)), collapse = " + "), 
    "+", random_effect
  ))
  
  # Fit the linear mixed-effects model
  model <- lmer(formula, data = papisQ_30)
  
  # Extract summary information
  model_summary <- summary(model)
  
  # Get coefficient, standard error, and p-value for the current variable
  coeff <- model_summary$coefficients[var, "Estimate"]
  std_err <- model_summary$coefficients[var, "Std. Error"]
  p_val <- model_summary$coefficients[var, "Pr(>|t|)"]
  
  # Calculations for combined effects at 10 years and 15 years
  # For 10 years
  sq10 <- sqrt(11)
  lg10 <- log(11)
  comb_10 <- paste(var, "+", sq10, "*", paste("s.time", ":", var), "+", lg10, "*", paste("l.time", ":", var), "= 0")
  mod_10 <- glht(model, linfct = comb_10)
  coef_10 <- summary(mod_10)$test$coefficients[[1]]
  std_err_10 <- summary(mod_10)$test$sigma[[1]]
  p_val_10 <- summary(mod_10)$test$pvalues[[1]]
  
  # For 15 years
  sq15 <- sqrt(16)
  lg15 <- log(16)
  comb_15 <- paste(var, "+", sq15, "*", paste("s.time", ":", var), "+", lg15, "*", paste("l.time", ":", var), "= 0")
  mod_15 <- glht(model, linfct = comb_15)
  coef_15 <- summary(mod_15)$test$coefficients[[1]]
  std_err_15 <- summary(mod_15)$test$sigma[[1]]
  p_val_15 <- summary(mod_15)$test$pvalues[[1]]
  
  # Return results as a data frame
  return(data.frame(
    variable = var,
    coefficient = coeff,
    std_error = std_err,
    p_value = p_val,
    coefficient_10 = coef_10,
    std_error_10 = std_err_10,
    p_value_10 = p_val_10,
    coefficient_15 = coef_15,
    std_error_15 = std_err_15,
    p_value_15 = p_val_15,
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
#                         Set up Parallel Computing
# =============================================================================

# Detect the number of cores and reserve some for other processes
num_cores <- detectCores() - 3  # Use all but one core

# Parallel computation using mclapply
results_list <- mclapply(loop_vars, fit_model, mc.cores = num_cores)

# Combine the results into a single data frame
results <- do.call(rbind, results_list)

# Preview the results and check the dimensions
head(results, 3)
dim(results)

# Save the results to a text file
write.table(results, paste(de, "/Direct_MLM_gwas_results_split_103.txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)


# =============================================================================
#                         GWAS Results Processing
# =============================================================================

# Clean the environment
rm(list = ls())

# Get the paths of the GWAS result files
paths_gwas <- fs::dir_ls("/gwas_results/")
paths_gwas <- paths_gwas[grepl("Direct", paths_gwas)]  # Filter files containing 'Direct'

# Check the number of GWAS result files
length(paths_gwas)

# Read the GWAS results from each file
data_gwas <- paths_gwas %>% map(function(path) {
  read.table(path, header = TRUE)
})

# Clean up file names for easier reference
names_gwas <- gsub("(.*/\\s*(.*$))", "\\2", paths_gwas)
names_gwas <- sub(".[^.]+$", "", names_gwas)

# Assign cleaned-up names to the data_gwas list
data_gwas <- set_names(data_gwas, names_gwas)

# Preview the data
head(data_gwas)


# =============================================================================
#                     Combining GWAS Data at SPC Level 1 Chronic Pain
# =============================================================================

# Combine all individual GWAS results into one data frame
data_gwas_comb <- as.data.frame(rbindlist(data_gwas))

# Check the dimensions and preview the combined data
dim(data_gwas_comb)
head(data_gwas_comb, 2)

# Sort the combined data by p-value
data_gwas_comb <- data_gwas_comb[order(data_gwas_comb$p_value),]

# Extract chromosome and position information from the SNP variable
data_gwas_comb$CHR <- as.numeric(gsub("chr", "", sub("_.*", "", data_gwas_comb$variable)))
data_gwas_comb$POS <- as.numeric(sub("_.*", "", sub("^[^_]*_", "", data_gwas_comb$variable)))
data_gwas_comb$SNP <- gsub("_", ":", sub("_[^_]*$", "", data_gwas_comb$variable))

# Preview the modified data
head(data_gwas_comb, 2)

# =============================================================================
#                              Plotting Manhattan Plot
# =============================================================================

# Load the qqman package for creating the Manhattan plot
library(qqman)

# Save the Manhattan plot as a PDF file
ggsave(
  plot = manhattan(
    data_gwas_comb,
    chr = "CHR",
    bp = "POS",
    p = "p_value_15",
    snp = "SNP",
    col = c("skyblue", "tomato"),  # Colors for different chromosomes
    suggestiveline = -log10(1e-5),  # Suggestive significance threshold
    genomewideline = -log10(5e-8),  # Genome-wide significance threshold
    main = "MLM GWAS",  # Title of the plot
    annotateTop = TRUE,  # Annotate the top SNPs
    ylim = c(0, 10),  # Set y-axis limits
    cex = 0.6  # Adjust text size
  ),
  filename = "manhattan_plotb.pdf",
  width = 10,
  height = 6
)
