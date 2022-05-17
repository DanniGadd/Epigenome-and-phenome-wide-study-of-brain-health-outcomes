################################################################################################################################

### SENSITIVITY ICV PREPS 

################################################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/

screen

R

library(tidyverse)
library(readxl)

# This script will format and combine all variables required for PheWAS analyses 
# These variables (N=1065 total) will be loaded into lmekin models for association tests
# Each variable will be tested for 4235 protein levels to identify protein markers of brain health outcomes

################################################################################################################################

### PROTEIN DATA

# Protein data (raw, without processing)
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/GS+soma+QC+normalized.csv")

# Annotation linker file for SeqIds
link <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/annotation.csv", check.names = F)

# Update naming so were working in SeqIds 
names <- colnames(link)
names(prot)[c(33:4267)] <- names

# Load target file 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/STRADL_DNAm_target_REM_17April2020.txt")

## TRANSFORM PROTEIN DATA AND JOIN FOR REGRESSIONS

## Log Transform 
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- log(prot[,i])
}


## Rank-Inverse Based Normaliation
library(bestNormalize)
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- orderNorm(prot[,i])$x.t
}

### SENSITIVITY CHECK FOR STUDY SITE AND LAG GROUP VARIABLES ON PROTEIN LEVELS:

# # Get lag group and study site variables 
# ewas <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/soma_demo1_file_for_rosie.csv")
# var <- ewas[c("st","study_site", "lag_group")]
# names(var)[1] <- "SampleId"

# # Add them into the main protein file 
# prot <- left_join(prot, var, by = "SampleId")

# # > dim(target)
# # [1] 847   6

# # Check for missing protein data 
# # proteins <- prot[c(13,33:4267)]

# ## TRANSFORM PROTEIN DATA AND JOIN FOR REGRESSIONS

# ## Log Transform 
# for(i in colnames(prot)[33:4267]){ 
#   prot[,i]<- log(prot[,i])
# }

# ## Regress Proteins onto Covariates 
# for(i in colnames(prot)[33:4267]){ 
#   prot[,i]<- lm(prot[,i] ~ factor(study_site) + factor(lag_group), 
#                                     na.action = na.exclude, data = prot)$residuals
# }

# ## Rank-Inverse Based Normaliation
# library(bestNormalize)
# for(i in colnames(prot)[33:4267]){ 
#   prot[,i]<- orderNorm(prot[,i])$x.t
# }

# prot_OG <- prot 
# prot_resid <- prot 

# table <- data.frame(SeqId = 1:4235, Corr = 1:4235)

# names <- as.character(colnames(prot)[33:4267])

# for (i in 1:length(names)){
#   marker <- names[i]
#   result <- cor.test(prot_OG[,marker], prot_resid[,marker])
#   estimate <- result$estimate
#   table[i,1] <- marker
#   table[i,2] <- estimate 
# }

# table <- table[order(table$Corr),]
# write.csv(table, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/table_correlations.csv", row.names = F)

## WMost proteins have a correlation of 95% or above, however a few are 0.87 and above
# For this reason, we will add study site and lag group into cognitive and apoe models as covariates 
# We will also add lag group into imaging models (as study_site*ICV is already adjusted for)

################################################################################################################################

### ADD PHENOTYPES

# Load demographics for all people in STRADL 
demo <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/demographicsV2.csv")
names(demo)[1] <- "SampleId"

# Join phenotype info to protein dataset 
prot <- left_join(prot, demo, by = "SampleId")

# Join in the GS id linkage so that further phenotypes can be joined in 
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/ST id linkage.csv")
names(IDs)[2] <- "GS_id"
names(IDs)[1] <- "SampleId"
IDs <- IDs[c(1,2)]
prot <- left_join(prot, IDs, by = "SampleId")

# Join in APOE
APOE <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/apoe_GS.csv")
names(APOE)[1] <- "GS_id"
prot <- left_join(prot, APOE, by = "GS_id")

# Check how many have each of the apoe phenotypes 
test <- prot
test$apoe %>% unique() #  "e3e4" "e3e3" "e2e3" NA     "e4e4" "e2e4" "e2e2"
table(is.na(test$apoe)) # 15 missing NA values 
outcomes <- test$apoe %>% as.data.frame() 
names(outcomes)[1] <- "X"
count(outcomes, X)

prot <- prot %>% mutate(APOE = case_when(
  apoe == "e4e4" ~ 2,
  apoe == "e3e4" ~ 2,
  apoe == "e3e3" ~ 1,
  apoe == "e2e2" ~ 0,
  apoe == "e2e3" ~ 0))

table(prot$APOE)
#   0   1   2
# 126 633 269

# Now join in the processed cognitive data from the script daniel shared with me - composite gf and g scores and other scores with outliers > 3.5 sd from mean removed 
comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/cog1_270121.csv")
names(comp)[1] <- "SampleId"
prot <- left_join(prot, comp, by = "SampleId")

# Add depression status into the dataset (read in file generated in depression covariate check, that is prepped with combined case/controls) 
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined_SCID_only_270121.csv")
dep <- dep[c(1,2)]
dep$DiagnosisGiven <- as.character(dep$DiagnosisGiven)
names(dep)[1] <- "SampleId"
prot <- left_join(prot, dep, by = "SampleId")
table(is.na(prot$DiagnosisGiven))
table(prot$DiagnosisGiven)

# Add BMI as a covariate 
demo_table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/protAA/demo_table.csv")
BMI <- demo_table[c(14,54)]
prot = left_join(prot, BMI, by = "GS_id")
library(imputeTS)
prot$BMI = na_mean(prot$BMI) 

# # Join in the GS id information into the protein dataset for models
# ids <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/ST id linkage.csv")
# names(ids)[1] <- "SampleId"
# ids <- ids[c(1,2)]
# prot <- left_join(prot, ids, by = "SampleId")


### ADD IN IMAGING DATA 

# Read in variables
brain_age <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_BrainAge.csv")
age <- brain_age[c(1,5)]
names(age)[1] <- "ID"

vol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_Volumes.csv")

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Updated_WM/STRADL_Brain_Measures_gFA-gMD.csv")

WMHV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_WMHV_110821/STRADL_WMHV_Measures_Complete_07.2021.csv")

faz <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_Fazekas.csv")

ICV_standardised <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_IDP_Measures_Standardised_ICV_Main.csv")

### Join
vol <- full_join(vol, WM, by = "ID", all = TRUE)
vol <- full_join(vol, WMHV, by = "ID", all = TRUE)
vol <- full_join(vol, age, by = "ID", all = TRUE)
vol <- full_join(vol, faz, by = "ID", all = TRUE)
vol <- full_join(vol, ICV_standardised, by = "ID", all = TRUE)
write.csv(vol, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/imaging_joint_ICV_added.csv", row.names = F)

# Recode points > 3.5 sd from mean for each variable (this has already been done for cognitive variables in the cognitive prep script for g and gf)
low = mean(vol$Global_GM_Volume, na.rm=T) - 3.5*sd(vol$Global_GM_Volume, na.rm=T)
high = mean(vol$Global_GM_Volume, na.rm=T) + 3.5*sd(vol$Global_GM_Volume, na.rm=T)
table(vol$Global_GM_Volume < low | vol$Global_GM_Volume > high)
vol$Global_GM_Volume[vol$Global_GM_Volume < low | vol$Global_GM_Volume > high] <- NA # no outliers GGM

low = mean(vol$WBV_No_Ventricles, na.rm=T) - 3.5*sd(vol$WBV_No_Ventricles, na.rm=T)
high = mean(vol$WBV_No_Ventricles, na.rm=T) + 3.5*sd(vol$WBV_No_Ventricles, na.rm=T)
table(vol$WBV_No_Ventricles < low | vol$WBV_No_Ventricles > high)
vol$WBV_No_Ventricles[vol$WBV_No_Ventricles < low | vol$WBV_No_Ventricles > high] <- NA # one outlier WBV

low = mean(vol$gFA, na.rm=T) - 3.5*sd(vol$gFA, na.rm=T)
high = mean(vol$gFA, na.rm=T) + 3.5*sd(vol$gFA, na.rm=T)
table(vol$gFA < low | vol$gFA > high)
vol$gFA[vol$gFA < low | vol$gFA > high] <- NA # one outlier gFA

low = mean(vol$gMD, na.rm=T) - 3.5*sd(vol$gMD, na.rm=T)
high = mean(vol$gMD, na.rm=T) + 3.5*sd(vol$gMD, na.rm=T)
table(vol$gMD < low | vol$gMD > high)
vol$gMD[vol$gMD < low | vol$gMD > high] <- NA # one outlier gMD

low = mean(vol$Brain_age, na.rm=T) - 3.5*sd(vol$Brain_age, na.rm=T)
high = mean(vol$Brain_age, na.rm=T) + 3.5*sd(vol$Brain_age, na.rm=T)
table(vol$Brain_age < low | vol$Brain_age > high)
vol$Brain_age[vol$Brain_age < low | vol$Brain_age > high] <- NA # no outliers brain age 

low = mean(vol$Fazekas_Score_Total, na.rm=T) - 3.5*sd(vol$Fazekas_Score_Total, na.rm=T)
high = mean(vol$Fazekas_Score_Total, na.rm=T) + 3.5*sd(vol$Fazekas_Score_Total, na.rm=T)
table(vol$Fazekas_Score_Total < low | vol$Fazekas_Score_Total > high)
vol$Fazekas_Score_Total[vol$Fazekas_Score_Total < low | vol$Fazekas_Score_Total > high] <- NA # 4 outliers

# +1 and log transform the WMHV variable
vol$WMH_Volume_Total <- as.numeric(vol$WMH_Volume_Total)
vol$WMH_Volume_Total <- vol$WMH_Volume_Total + 1
vol$WMH_Volume_Total <- log(vol$WMH_Volume_Total)

low = mean(vol$WMH_Volume_Total, na.rm=T) - 3.5*sd(vol$WMH_Volume_Total, na.rm=T)
high = mean(vol$WMH_Volume_Total, na.rm=T) + 3.5*sd(vol$WMH_Volume_Total, na.rm=T)
table(vol$WMH_Volume_Total < low | vol$WMH_Volume_Total > high)
# vol$gMD[vol$WMH_Volume_Total < low | vol$WMH_Volume_Total > high] <- NA # 6 outlying values for WMHV
vol$WMH_Volume_Total[vol$WMH_Volume_Total < low | vol$WMH_Volume_Total > high] <- NA # 6 outlying values for WMHV

### Join the imaging data into the protein dataset for regressions 
# We want imaging data for all of the individuals with protein data, so will leftjoin 
names(vol)[1] <- "SampleId"
prot <- left_join(prot, vol, by = "SampleId")

# Calculate brain age acceleration score 
prot$brain_accel <- resid(lm(Brain_age ~ st_age, na.action = na.exclude, data = prot))

# add in lag group and study site 
ewas <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/soma_demo1_file_for_rosie.csv")
var <- ewas[c("st","study_site", "lag_group")]
names(var)[1] <- "SampleId"

# Add them into the main protein file 
prot <- left_join(prot, var, by = "SampleId")

# # Save out prot file that is prepped wth proteins and phenotypes
# write.csv(prot, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", row.names = F)

###################

### LOAD COX REQUIREMENTS 


## Load lmekin requirements
library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

## Code that is already processed and read in as the ped file below
# ped = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree.csv")
# ped$father <- as.numeric(ped$father)
# ped$mother <- as.numeric(ped$mother)
# ped$father[ped$father==0] <- NA
# ped$mother[ped$mother==0] <- NA
# table(ped$sex)
# ped$sex <- as.numeric(ped$sex)
# ped$sex[ped$sex==2] <- 0
# ped$sex <- ped$sex+1

# Read in the prepped file to cluster 
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")

# Create kinship matrix for GS
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid)) # Pedigree list with 26 total subjects in 5 families
kin_model <- kinship(kin) 

# Function to Extract Lmekin Results  
extract_coxme_table <- function (mod){
  #beta <- mod$coefficients #$fixed is not needed
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- beta/se
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

## Set marker names to loop through in models as proteins (x)
markers <- prot[c(33:4267)]
# marker_names <- colnames(markers)

# Rename depression variable
names(prot)[4297] <- "combined"

####################################

### SENSITIVITIES: GLOBAL GREY MATTER VOLUME 

# Set marker names to top assocs from study 
GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Global_GM_Volume_result_annotated_THR.csv")
GGM$type <- "imaging"
GGM$phenotype <- "Global Grey Matter Volume"
GGM <- GGM[which(GGM$Status == "pass"),]
marker_names <- GGM$SeqId

markers <- markers[,which(colnames(markers) %in% marker_names)]

# Generate results file with the standard approach we have used 
# i.e. ICV*Stduy_site without precorrecting ICV by study site first, then + batch + editor as volume specific covariates 
# The sig assocs are (STC1, SLITRK1, IGLON5, ANG, NCAN, NEU1, PTPRD, CFHR1, TREM1, RBL2, C5, SEZ6L, NPTXR, TNFRSF11B, TPP1)

length <- 15
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:15){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/ICV_sensitivity/", pheno, "_result_as_per_paper.csv"))
}



# Generate results file for GGM/corrected_ICV + study_site + batch + editor in models 

prot$GGM_adjusted <- prot$Global_GM_Volume / prot$Standardised_ICV

length <- 15
phenotype <- c("GGM_adjusted")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:15){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/ICV_sensitivity/", pheno, "_result_GGM_over_ICV_corrected_site_batch_edited.csv"))
}



# Generate results file for GGM/ICV + study_site + batch + editor in models 


prot$GGM_adjusted_2 <- prot$Global_GM_Volume / prot$Estimated_ICV

length <- 15
phenotype <- c("GGM_adjusted_2")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:15){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/ICV_sensitivity/", pheno, "_result_GGM_over_ICV_uncorrected_site_batch_edited.csv"))
}


# Generate results for corrected_ICV*study_site _ batch + editor in models 

length <- 15
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:15){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Standardised_ICV + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/ICV_sensitivity/", pheno, "_result_as_per_paper_but_with_standardised_ICV.csv"))
}


# Generate results as per study but without the interaction term and with standardised ICV (i.e. all separate covariates)

length <- 15
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:15){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + prot$Standardised_ICV + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/ICV_sensitivity/", pheno, "_result_as_per_paper_but_without_interaction_term_and_with_standardised_ICV.csv"))
}


# Generate results as per study but without the interaction term and without standaridising ICV (i.e. all separate covariates)

length <- 15
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:15){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + prot$Estimated_ICV + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/ICV_sensitivity/", pheno, "_result_as_per_paper_but_without_interaction_term_and_without_standardised_ICV.csv"))
}



# Take a look at results 
# It seems like using standardised vs unstandardised ICV has not made much difference when these are included as covariates (either as interaction term or singularly)
# However, dividing GGM volume by the ICV variables first before running models seems to have a massive effect (i.e. many of the previously significant 15 assocs no longer significant)








