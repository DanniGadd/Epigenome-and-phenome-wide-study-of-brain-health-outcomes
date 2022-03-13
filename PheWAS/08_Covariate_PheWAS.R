################################################################################################################################

### PheWAS analyses - Smoking, BMI

################################################################################################################################

cd 

screen

R

library(tidyverse)
library(readxl)

## Load prepped joint PheWAS phenotypes file 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv", check.names = F)

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
marker_names <- colnames(markers)

# Rename depression variable
names(prot)[4297] <- "combined"

#####################################################################################################################################

### BMI

#####################################################################################################################################

length <- 4235
phenotype <- c("BMI")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/06_Covariates/batches/", pheno, "_result.csv"))
}

#####################################################################################################################################

### SMOKING

#####################################################################################################################################

length <- 4235
set <- c("set")
phenotype <- c("EverSmoker")

for (j in 1:length(set)){
  pheno <- "EverSmoker"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot[,pheno]) + as.factor(prot$combined) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
   write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/06_Covariates/batches/", pheno, "_result.csv"))
}






