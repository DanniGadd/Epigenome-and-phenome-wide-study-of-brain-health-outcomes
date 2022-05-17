############################################################################################

### MOA 774 with genetic ORM

############################################################################################

# Run 31/01/21

# EWAS with MOA on 774 individuals, using a genetic ORM to adjust for relatedness
# Adjustement for basic covariates:
# age
# sex
# wave
# batch 
# Depression status (enriched in STRADL)
# EXTRA: the pQTLs have been regressed from proteins

# These models will produce the pQTL-adjusted basic results.

# The DNAm file from the basic script willbe used for this script, as the DNAm does not change (only 
# the pQTLs added to the protein regression)

############################################################################################################

### Phenotype file - 778

############################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/01_Basic

screen

R

# Read in order of meth files we are matching to - combined IDS of individuals (from the relatedness subtyping script for GRM)
order <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")

# Read in protein file which has GS id and stradl id, then all 4,235 proteins that have been preprocessed
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Protein_preps/Phenotype_file_774_pQTLs_regressed.csv", check.names = F)
names(pheno2)[1] <- "ID"
pheno2 <- pheno2[-2]

# Check to make sure the id in the proteins file is the same as the id in the ID reference file 
identical(pheno2$ID, order$GS_id) # TRUE

# Assign FID and IID from the ID reference file 
pheno2$IID <- order$IID
pheno2$FID <- order$FID

# Order so proteins are first 
pheno2 <- pheno2[c(2:4236,1,4237,4238)]

identical(pheno2$IID, order$IID) # TRUE

# Write out phenotype files for a heritability run in the 4235
location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/"

## Write out protein files so FID, IID, phenotype 
for(i in 1:350){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch1/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 351:700){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch2/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 701:1050){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch3/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1051:1400){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch4/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1401:1750){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch5/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1751:2100){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch6/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2101:2450){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch7/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2451:2800){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch8/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2801:3150){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch9/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3151:3500){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch10/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3501:3850){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch11/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3851:4235){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch12/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

# The phenotype (protein) files are prepped per protein in the 12 batch folders to run in parallel

############################################################################################################

### BINARY FILES AND ORM - 774

############################################################################################################

# Use the same DNAm binary files and opi files as generated in the basic model script without pQTL adjustment
# This is because only the proteins have changed (been regressed for pQTLs)

##########################################################################################

# BATCH EWAS 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch1/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch2/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch3/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch4/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch5/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch6/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch7/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch8/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch9/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch10/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch11/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch12/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_basic_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 



####################################################################################################

### PROCESSING THE EWAS RESULTS FOR BATCHES ABOVE - MOA 778 with genetic ORM

####################################################################################################

# Split this up for speed into 12 screens 

# Go through each protein and subset to the significant sites 
# First i do the results cut to the standard EWAS significance threshold (without correction for multiple testing)

screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch1/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b1")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein 
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch2/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b2")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch3/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b3")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch4/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b4")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}




screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch5/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b5")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch6/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b6")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch7/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b7")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch8/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b8")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch9/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b9")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch10/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b10")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch11/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b11")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}


screen

R


library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/batch12/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

batch <- c("b12")

for(i in 1:length(list)){

	path <- list[i]
	batch <- batch[i]
	setwd(paste0(path))

	L <- list.files(".", ".mlma")

	# results tables for each protein (all vs sig)
	for(j in 1:length(L)){
		file <- fread(L[j], header = T)
		file <- as.data.frame(file)
		name <- as.character(L[j])
		name <- gsub("_.*", "", name)
		file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
		write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
	}
}

##########################################################################

# Collate top hits for proteins 

screen

R

path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

setwd(path)

L <- list.files(".", ".csv")
L # 4231 converged 

files <- lapply(L, read.csv)
names <- as.character(L)
batch <- gsub("_.*", "", names)
marker <- gsub("_results.*", "", names)
marker <- gsub(".*_", "", marker)
names(files) <- marker
osca <- do.call(rbind, files)
osca <- osca[c(9,1:8)]

# annotations file add in
library(tidyverse) 
library(readxl)
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")
comb <- left_join(osca, anno, by = "SeqId")


# Filter to CpG-protein correction value of adjusted p 0.05/143/772619 = 4.525521e-10
# 0.0000000004525521


comb2_filt2 <- comb[which(comb$p < 0.0000000004525521),]

dim(comb2_filt2) # 238245  - at cpg/protein adjusted significance 
comb2 <- comb2_filt2[order(comb2_filt2$p),]
write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_second_threshold_020221.csv", row.names = F)


# Filter to CpG-protein correction value of adjusted p 0.05/4235/772619 = 1.528098e-11
# > 0.05 / 4235 / 772619
# [1] 1.528098e-11

# 0.00000000001528098

comb2_filt2 <- comb[which(comb$p < 0.00000000001528098),]

dim(comb2_filt2) # 158858  - at cpg/protein adjusted significance 
comb2 <- comb2_filt2[order(comb2_filt2$p),]
write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_020221.csv", row.names = F)

##############################################################

### WORK OUT WHICH FULL PROTEINS DIDNT CONVERGE 

# Read in protein file which has GS id and stradl id, then all 4,235 proteins that have been preprocessed
pheno1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_4235_778_eGFR_included.csv", check.names = F)
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_removed.csv", check.names = F)
names(pheno2)[1] <- "ID"
pheno2 <- pheno2[-2]
colnames(pheno2) <- colnames(pheno1)


path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/02_Basic_pQTLs/batches/sig_hits/"

setwd(path)

L <- list.files(".", ".csv")
L

files <- lapply(L, read.csv)
names <- as.character(L)
batch <- gsub("_.*", "", names)
marker <- gsub("_results.*", "", names)
marker <- gsub(".*_", "", marker)
names(files) <- marker
osca <- do.call(rbind, files)
osca <- osca[c(9,1:8)]

# load pheno2 above and do this to work out difference 
list <- names(files)
list2 <- pheno2[-which(colnames(pheno2) %in% list)]

# Plots - of proteins look okay 
proteins <- colnames(list2)[1:8]
proteins


# NON CONVERGANCE pQTLs:

#     ID    15509-2    15584-9  15585-304   15587-20    15588-17     15589-1
# 1  243 -0.7622209 -1.7307952 -0.4091063 -0.4021089 -0.09844514  0.37776971
# 2  414  1.4475706  0.5464439  0.5464439  0.8561045 -0.66654991 -0.93302452
# 3  520  1.1244222 -1.0377357 -0.2357978 -0.1049236 -0.82400649  0.09197085
# 4  857  1.0432717 -1.2988312 -0.5167521 -1.2551504 -0.25572573  0.96332820
# 5 2199  1.1679172 -1.2411103  1.6758513 -0.2324858 -0.78394278  0.38817257
# 6 2238  1.6255440  0.7927355  0.6746192  0.9480676 -0.36740756  1.22049757
#      15591-28    15594-47     15596-7   15602-43   15603-20    15604-18
# 1  0.44796340 -0.87013335  1.21374068 -1.9171161 -0.7752102  0.94302952
# 2 -0.14064109  0.58425234 -0.49114108 -0.9633282  0.4266886 -0.24907219
# 3  0.28245661  1.10049555 -0.80604083 -0.2324858 -1.7757188 -0.39164948
# 4  0.08550041  0.04674334  0.01127901  0.3777697  0.7752102 -0.06288025
# 5 -1.27668647  0.02094782 -1.77571881 -1.1305052 -1.8972680 -1.33722196
# 6 -1.82453980  0.92311201 -0.02739472  0.8891094 -0.5919135 -0.47665058
#      15606-19   15607-56    15608-5     15612-5    15613-16  15614-168
# 1 -0.27575570 -1.0771826  0.3639623  0.04351755  0.17981298 -1.3372220
# 2 -1.42054977 -1.0714452 -0.6505399 -0.98927368 -0.16019632 -0.2391124
# 3 -0.29253190 -0.5167521 -0.7752102 -1.16156870  0.08226655  0.1504115
# 4 -0.06288025  0.6386410  0.2457496 -0.23248578  0.69909604 -1.2273106
# 5 -0.39861759  0.9280569 -1.3945277 -1.63773578  0.85610447  0.1471532
# 6 -0.76222094 -1.1872509 -0.4408497 -0.21596306  1.54628259  0.1049236
#      15615-8     15617-8     15620-4    15622-13    4407-10     6402-8
# 1  0.3127730  0.03061858 -1.03223127 -1.13050522  2.1121964 -0.2258693
# 2 -0.3467997 -0.68273261  0.30939081  0.02094782  0.7239989  0.2925319
# 3 -0.7622209 -1.41176934 -0.21596306  0.40560510  0.3847003  0.8561045
# 4  0.2657274 -1.52533171 -0.91818962  0.43730125  0.6307584 -0.5539417
# 5 -0.4730437 -1.61358906 -0.74085280  0.83308965 -1.4205498 -0.2757557
# 6 -0.3674076 -0.30601214  0.04029221  0.32294129 -0.3229413  0.2724099

