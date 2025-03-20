#_______________________________________________________________________________
#        Multilevel model for genotypes
#_______________________________________________________________________________
library(tidyr);library(lattice); library(lmtest); library(readxl);library(car); library(backports)     # to revive the isFALSE() function for sim_slopes()
library(effects);library(ggplot2);library(interactions); library(lme4); library(lmerTest)      
library(psych); library(plyr);library(dplyr);library(R2MLwiN);library(data.table)
options(MLwiN_path='/Users/uzoemeka/Documents/mlnscript')
library(writexl);library(haven);library(tidyverse);library(nlme);library(curl);library(rstan);library(brms)
library(lubridate);library(kableExtra);library(foreign);library(multcomp);library(groupdata2);library(merTools)


#====================================================================================
#                 DATA
#====================================================================================
rm(list = ls())
papis <- read_dta("/Users/uzoemeka/Documents/MLM_Analysis/MLM/Edssdatatata2nd_three_month.dta")
papis <- as.data.frame(papis)
names(papis)[names(papis) == 'date'] <- 'relapseDate'

#______________________________________________________________________________
papis$patid <- papis$patientIdentifier
papis <- papis[,c("patid","dateTime","msDateDiagnosis", "msDateOnset", "DMT",                                     
                  "relapseDate", "edssScore", "edssScoreDecimal", "msEventType",
                  "sex", "dob")] 


#______ Make R recognise other date variable using as.Date                                                            
papis$msDateDiagnosis <- as.Date(papis$msDateDiagnosis, "%d/%m/%Y")                                         
papis$msDateOnset <- as.Date(papis$msDateOnset, "%d/%m/%Y")
papis$dateTime <- as.Date(papis$dateTime, "%d/%m/%Y")
papis$dob <- as.Date(papis$dob, "%d/%m/%Y")
papis$date <- as.Date(papis$relapseDate, "%d/%m/%Y")

#_______________________ Create timeto and age _________________
papis$ageonset <- as.numeric((papis$msDateOnset - papis$dob)/365.25)
papis$timeage <- as.numeric((papis$dateTime - papis$dob)/365.25)
papis$timediag <- as.numeric((papis$dateTime - papis$msDateDiagnosis)/365.25)
papis$timeonset <- as.numeric((papis$dateTime - papis$msDateOnset)/365.25)

#papis[papis$timeonset < 0, ] 
papis <- papis[papis$timeonset >= 0,]
#papis[papis$timeonset < 0, ] 

papis <- papis %>% mutate (year_diagnosis=as.numeric(format(msDateDiagnosis, '%Y')))
papis$year_diagnosis_cat <- as.factor(ifelse(papis$msDateDiagnosis > as.POSIXct("1953", format = "%Y") & papis$msDateDiagnosis < as.POSIXct("2005", format = "%Y") , "year1",
                                             ifelse(papis$msDateDiagnosis > as.POSIXct("2006", format = "%Y") & papis$msDateDiagnosis < as.POSIXct("2010", format = "%Y"), "year2", "year3")))

#_____________________ Create nonlinear time terms _________________
papis$time <- papis$timeonset + 1
papis$l.time <- log(papis$time)
papis$s.time <- sqrt(papis$time)

papis$DMT <- as.factor(papis$DMT)
papis$DMT <- relevel(papis$DMT, ref = "notDMT")

papis$sex <- as.factor(papis$sex)
papis$sex <- relevel(papis$sex, ref = "MALE")

#______________ Create categorical sex and DMT and year of diagnosis_________________
papis$sex.cat <- ifelse(papis$sex == "MALE", 0, 1)
papis$DMT.cat <- ifelse(papis$DMT == "notDMT", 0, 1)

papis$year.diagnosis.cat <- ifelse(papis$year_diagnosis_cat == "year1", 0,
                                   ifelse(papis$year_diagnosis_cat == "year2", 1, 2))

#_______________ Observation within 1/4 year period __________
id <- unique(papis$patid)
papis_quarter <- NULL

for(i in id){
  papis_i <- as.data.frame(papis %>% filter(patid == i) %>%
                             dplyr::select("patid", "dateTime", "edssScore",
                                           "ageonset", "timeage", "timeonset", "time",
                                           "s.time",  "l.time", "sex.cat", "DMT.cat", "year.diagnosis.cat")
  )
  break_dates <- seq.Date(min(papis$dateTime),
                          as.Date("2022-01-1"),
                          by = "quarter")
  cut_i <- cut(papis_i$dateTime, breaks = break_dates)
  papis_median <- aggregate(papis_i,by=list(cut_i),median)
  
  papis_quarter <- rbind(papis_quarter, papis_median)
}
papis_quarter$DMT.cat[papis_quarter$DMT.cat == 0.5]<- 1
length(unique(papis_quarter$patid))

pca <- fread("/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/karim_clean_threemonths_pca1.eigenvec")
pca <- pca[,c(1:7)]
pca$ids <- sub("^[^_]*_", "", pca$FID)
head(pca,2)
length(unique(pca$FID))


# No_year_and_dmt_top_10_15yr_imputed.raw
# karim_clean_threemonths_imputed_noFID_1to100k.raw
# karim_clean_threemonths_imputed_noFID_1to1k.raw
# karim_clean_threemonths_imputed_noFID_500kto5001k.raw


di <- "/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/Raw_data_split_3"
raw_file <- fread(paste(di, "/karim_clean_threemonths_imputed_noFID_split_103.raw", sep = ""))
head(raw_file[,1:10],2); dim(raw_file)
#raw_file <- raw_file[,1:10]
raw_file$ids <- sub("^[^_]*_", "",sub("^[^_]*_", "",sub("^[^_]*_", "",sub("^[^_]*_", "", raw_file$IID))))
length(unique(raw_file$IID))
length(unique(raw_file$ids))

ids <- fread("/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/all.ids_Genotyped_patients_threemonthspostrelapse.csv")
ids$ids <- ids$IID

rawfile <- merge(raw_file, ids, by = "ids")
length(unique(rawfile$ids))
# names(rawfile)
# rawfile <- rawfile[,-c( 'IID.y', 'IID.x' , 'sex.cat', 'ageonset', 'ageonset.cent', 'dupl.x', 
#                         'FID', 'PatID', 'MatID', 'SEX', 'PHENO', 'dupl.y')]

# rawfile <- rawfile[,-c("IID.x", "PAT", "MAT", "SEX.x","PHENOTYPE","IID.y",
#                         "sex.cat", "ageonset","ageonset.cent", "dupl.x","FID", "PatID", "MatID",  "SEX.y",
#                         "PHENO","dupl.y")]

rawfile <- rawfile[,-c("IID.x","IID.y",
                       "sex.cat", "ageonset","ageonset.cent", "dupl.x","FID", "PatID", "MatID",  "SEX",
                       "PHENO","dupl.y")]


rawfile <- rawfile %>% mutate(across(.cols = -all_of(c("patid","ids")), .fns = scale))
colnames(rawfile)[-c(1,dim(rawfile)[2])] <- gsub(":", "_", colnames(rawfile)[-c(1,dim(rawfile)[2])])
head(rawfile[,1:10],3)
dim(rawfile)[2]
length(unique(rawfile$ids))

rawfile1 <- merge(rawfile, pca, by = "ids")
rawfile1 <- rawfile1[,-c( "ids", 'IID', 'FID')]
head(rawfile1[,1:10],2); dim(rawfile1)
length(unique(rawfile1$patid))

papis_quarter <- merge(papis_quarter, rawfile1, by = "patid")
head(papis_quarter[,1:10],2);dim(papis_quarter)
length(unique(papis_quarter$patid))



#____________________ Defining other covariates __________________________________________
papis_quarter$sex.cat <- factor(papis_quarter$sex.cat, labels = c("male", "female"))
papis_quarter$DMT.cat <- factor(papis_quarter$DMT.cat, labels = c("offDMT", "onDMT"))
papis_quarter$year.diagnosis.cat <- factor(papis_quarter$year.diagnosis.cat, labels = c("year1", "year2", "year3"))
papis_quarter$yeardiagnosiscat <- factor(papis_quarter$year.diagnosis.cat, labels = c("year1", "year2", "year3"))

papis_quarter$ageonset.cent <- papis_quarter$ageonset - 31.5
final_data <- papis_quarter#[complete.cases(papis_quarter),]
length(unique(final_data$patid))


#_______________________________________________________________________________
#                     30 year truncated data
#_______________________________________________________________________________
papisQ_30 <- as.data.frame(final_data %>% group_by(patid) %>% filter(time <= 30))
length(unique(papisQ_30$patid))

papisQ_30$sexcat <- papisQ_30$sex.cat 
papisQ_30$DMTcat <- papisQ_30$DMT.cat
papisQ_30$stime <- papisQ_30$s.time
papisQ_30$ltime <- papisQ_30$l.time
papisQ_30$ageonsetcent <- papisQ_30$ageonset.cent
papisQ_30$yeardiagnosiscat <-papisQ_30$year.diagnosis.cat


#
# names(papisQ_30)
# Exclude mpg and am for example
variables <- names(papisQ_30)[!names(papisQ_30) %in% c("Group.1","ageonset","timeage","timeonset","time",
                                                       "s.time","l.time","dateTime","sex.cat","patid", "edssScore",
                                                       "ageonset.cent",  "year.diagnosis.cat", "DMT.cat","yeardiagnosiscat","sexcat",            
                                                       "DMTcat","stime", "ltime","ageonsetcent", "PC1","PC2","PC3","PC4","PC5")]  

# Load necessary libraries
library(parallel)  # For parallel computing
library(lme4)
library(broom.mixed)

# Define your data and model variables
constant_vars <- c("1", "s.time*ageonset.cent", "l.time*ageonset.cent", "s.time*sex.cat", "l.time*sex.cat", "PC1", "PC2", "PC3", "PC4", "PC5")
loop_vars <- variables  # Your looped variables
random_effect <- "(1 + s.time + l.time | patid)"  # Random effect variable

# Function to fit the model and extract the results for a single variable
fit_model <- function(var) {
  # Create the formula
  formula <- as.formula(paste("edssScore ~", paste(c(constant_vars, paste("l.time", "*", var), paste("s.time", "*", var)), collapse = " + "), "+", random_effect))
  
  # Fit the linear mixed effect model
  model <- lmer(formula, data = papisQ_30)
  
  # Extract coefficients, standard errors, and p-values for the current variable
  model_summary <- summary(model)
  coeff <- model_summary$coefficients[var, "Estimate"]
  std_err <- model_summary$coefficients[var, "Std. Error"]
  p_val <- model_summary$coefficients[var, "Pr(>|t|)"]
  
  # Calculations for 10 and 15-year combined effects
  sq10 <- sqrt(11)
  lg10 <- log(11)
  comb_10 <- paste(var, "+", sq10, "*", paste("s.time", ":", var), "+", lg10, "*", paste("l.time", ":", var), "= 0")
  mod_10 <- glht(model, linfct = comb_10)
  coef_10 <- summary(mod_10)$test$coefficients[[1]]
  std_err_10 <- summary(mod_10)$test$sigma[[1]]
  p_val_10 <- summary(mod_10)$test$pvalues[[1]]
  
  sq15 <- sqrt(16)
  lg15 <- log(16)
  comb_15 <- paste(var, "+", sq15, "*", paste("s.time", ":", var), "+", lg15, "*", paste("l.time", ":", var), "= 0")
  mod_15 <- glht(model, linfct = comb_15)
  coef_15 <- summary(mod_15)$test$coefficients[[1]]
  std_err_15 <- summary(mod_15)$test$sigma[[1]]
  p_val_15 <- summary(mod_15)$test$pvalues[[1]]
  
  # Return results as a list
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

# Set up parallel computing
num_cores <- detectCores() - 3  # Use all but one core

# Use mclapply for parallel execution
results_list <- mclapply(loop_vars, fit_model, mc.cores = num_cores)

# Combine the results into a single dataframe
results <- do.call(rbind, results_list)

head(results,3); dim(results)
de <- "/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/gwas_results"
write.table(results, paste(de, "/Direct_MLM_gwas_results_split_103.txt", sep = ""), col.names = T, row.names = F, quote = F)




#____________ GWAS results
rm(list=ls())
paths_gwas <- fs::dir_ls("/Users/uzoemeka/Documents/MLM_Analysis/MLM/MS_GWAS/Karim/gwas_results/")
paths_gwas <- paths_gwas[which(grepl("Direct",paths_gwas))]
length(paths_gwas)

data_gwas <- paths_gwas %>% map(function(path){
  read.table(path, header=TRUE)
})

names_gwas <- gsub("(.*/\\s*(.*$))", "\\2", paths_gwas)
names_gwas <- sub(".[^.]+$", "", names_gwas)


data_gwas <- set_names(data_gwas, names_gwas);head(data_gwas)


#________ SPC level 1 chronic pain combined 
data_gwas_comb <- as.data.frame(rbindlist(data_gwas))
dim(data_gwas_comb)
head(data_gwas_comb,2)
data_gwas_comb[order(data_gwas_comb$p_value),]

data_gwas_comb$CHR <- as.numeric(gsub("chr", "", sub("_.*", "", data_gwas_comb$variable)))
data_gwas_comb$POS <- as.numeric(sub("_.*", "", sub("^[^_]*_", "", data_gwas_comb$variable)))
data_gwas_comb$SNP <- gsub("_", ":", sub("_[^_]*$", "", data_gwas_comb$variable))
head(data_gwas_comb,2)
# str(data_gwas_comb)


library(qqman)
ggsave(plot = manhattan(data_gwas_comb,
                        chr = "CHR",
                        bp = "POS",
                        p = "p_value_15",
                        snp = "SNP",
                        col = c("skyblue", "tomato"),  # Colors for different chromosomes
                        suggestiveline = -log10(1e-5),  # Suggestive significance threshold
                        genomewideline = -log10(5e-8),  # Genome-wide significance threshold
                        main = "MLM GWAS",
                        annotateTop = TRUE,  # Annotate the top SNPs
                        ylim = c(0, 10),  # Setting y-axis limits
                        cex = 0.6  # Adjust text size
),
"manhattan_plotb.pdf", width = 10, height = 6)
