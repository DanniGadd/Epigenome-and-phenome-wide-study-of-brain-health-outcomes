
# So I've done the following, though I;m not sure how correct it is
# I managed to get an error-free output from it though
# 1 - Recode SNP chip data to 0,1,2 for the 778 individuals (IDs without GS suffix)
# 2 - Generate myprofile.txt from recoded dta
# 3 - Generate orm from -efile (myprofile.txt) using --make-orm
# 4 - Edit orm.id file in R, appending GS_ suffixes
# 5 - Run moa with DNAm data, fitting orm
# One difference between mine and yours is I just used the chip data - did you use the imputed data for your grm?
# I'm not so sure I've even done it right. OSCA documentation only really talks about gene expression/DNAm data
# I just substituted it with a matrix of 0,1,2 genotype counts
# They reference GWAS data for things like eQTL mapping

setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/")
# 758 samps
o778 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order.csv")

# 658 samps 
o658 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_658_meth_order.csv")

keep_778  = o778  
keep_778[,2] = gsub("GS_", "", o778[,2])
keep_778[,3] = gsub("GS_", "", o778[,3])
write.table(keep_778[,2:3], file="keep_778.txt", sep=' ', quote=F, row.names=F, col.names=F)

keep_658  = o658
keep_658[,2] = gsub("GS_", "", o658[,2])
keep_658[,3] = gsub("GS_", "", o658[,3])
write.table(keep_658[,2:3], file="keep_658.txt", sep=' ', quote=F, row.names=F, col.names=F)

# 778 SNP prep to orm

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/

# plink19 --recodeA --bfile /GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015 --out osca_test_778 --keep keep_778.txt

snps = fread("osca_test_778.raw", header=T, stringsAsFactors=F)
snps = as.data.frame(snps)
snps[,3:6] = NULL

write.table(snps, file="myprofile.txt", row.names=F, sep=' ')

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/
osca_Linux --efile myprofile.txt --keep keep_778.txt --make-orm --out ormtest778

# osca_Linux --moa --efile myprofile.txt --pheno $i --qcovar /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/WBC_covariates_778_main_model.cov --fast-linear --out test1


# 658 SNP prep to orm

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/
# plink19 --recodeA --bfile /GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015 --out osca_test_658 --keep keep_658.txt

snps = fread("osca_test_658.raw", header=T, stringsAsFactors=F)
snps = as.data.frame(snps)
snps[,3:6] = NULL

write.table(snps, file="myprofile658.txt", row.names=F, sep=' ')

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/
osca_Linux --efile myprofile658.txt --keep keep_658.txt --make-orm --out ormtest658

# osca_Linux --moa --efile myprofile658.txt --pheno $i --qcovar /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/WBC_covariates_658_main_model.cov --fast-linear --out test_658_daniel


# # Set pheno 

# i=/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_120421/test_1_778/2780-35eGFR_included_778.phen

# osca_Linux --moa --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order --pheno $i --qcovar /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/WBC_covariates_778_main_model.cov --fast-linear --out test_MOA_778_output --methylation-m

# 778 - recode the GS ids to match format in other files 

# In R (set IDs to GS_XXX to match DNAm IDs)
orm = read.table("ormtest778.orm.id")
orm[,1] = paste0("GS_", orm[,1])
orm[,2] = paste0("GS_", orm[,2])

write.table(orm, file="ormtest778.orm.id", sep='\t', col.names=F, row.names=F, quote=F)

# 658 - recode the GS ids to match format in other files 

# In R (set IDs to GS_XXX to match DNAm IDs)
orm = read.table("ormtest658.orm.id")
orm[,1] = paste0("GS_", orm[,1])
orm[,2] = paste0("GS_", orm[,2])

write.table(orm, file="ormtest658.orm.id", sep='\t', col.names=F, row.names=F, quote=F)
