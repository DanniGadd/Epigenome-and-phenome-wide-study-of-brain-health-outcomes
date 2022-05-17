
########################################################################

### READ IN NEURO CPGS AND PLOT ASSOCIATIONS WITH PROTEINS

screen 

R

slice <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_41_pQTMs_pQTLs_added.csv")

prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv", check.names = F)

# # subset methylation to cpgs of interest
# meth_neuro <- full[,which(colnames(full) %in% slice$CpG)]
# meth_neuro <- as.data.frame(meth_neuro)
# meth_neuro$GS_id <- rownames(meth_neuro)

# Join protei data to methylation data 
prot$GS_id <- as.character(prot$GS_id)
# neuro <- left_join(meth_neuro, prot, by = "GS_id")
# write.csv(neuro, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/plots_neuro_assocs/plots/neuro_joint_protein_cpgs.csv")

neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/plots_neuro_assocs/plots/neuro_joint_protein_cpgs.csv", check.names = F)

# Plot associations of interest 
table <- slice[c("CpG", "SeqId", "Gene.of.Protein")]

library(ggpubr)
library(ggplot2)

plots <- list()

for(i in 1:length(table$CpG)){
    CpG <- as.character(table[i,1])
    protein <- as.character(table[i,2])
    gene <- as.character(table[i,3])
    data <- neuro[,CpG] %>% as.data.frame()
    data2 <- neuro[,protein] %>% as.data.frame()
    dataset <- cbind(data,data2)
    names(dataset) <- c("CpG", "Protein")
    # pdf(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/plots_neuro_assocs/plots/", CpG, "_", protein, ".pdf"))
    #  ggplot(dataset, aes(x=CpG,y=Protein)) +
    #  geom_point(alpha=0.5) +
    #  labs(x= CpG, y=protein)
    # dev.off()
     p <- ggplot(dataset, aes(x=CpG,y=Protein)) +
     geom_point(alpha=0.5) +
     labs(x= CpG, y=gene)
    plots[[i]] <- p
}

# Save off
pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/plots_neuro_assocs/plots/neuro_pQTMs.pdf"))
for (i in 1:41) {
    print(plots[[i]])
}
dev.off()

############################################################################

### CREATE SUPPLEMENTARY FIGURES WITH STITCHED PLOTS FOR THE 41 pQTMs
# Split by cis/trans and group by protein marker 

# Plot associations of interest 
table <- slice[c("CpG", "SeqId", "Gene.of.Protein", "Effect")]

