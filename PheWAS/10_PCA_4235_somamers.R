##############################################################################################################

### PROCESS PheWAS models

##############################################################################################################

## Assess the independent signals in the 4235 somamers for inclusion in PheWAS/MWAS
## This will be used to decide the threshold for multuple testing adjustment in PheWAS vs FDR levels

# prcomp: https://www.analyticsvidhya.com/blog/2016/03/pca-practical-guide-principal-component-analysis-python/
# prcomp: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

screen

R

library(tidyverse)
library(ggplot2)
library(readxl)
library(psych)
library(ggcorrplot)
library(cowplot)
library(tidyverse)

# Read in the protein file 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv", check.names = F)

# Isolate just the protein columns of interest for PCA
prot <- prot[c(33:4267)]

# Try prcomp
library(factoextra)
res.pca <- prcomp(prot, scale = TRUE)

names(res.pca)
# [1] "sdev"     "rotation" "center"   "scale"    "x"

loadings <- res.pca$rotation # this provides the loadings

loadings[1:5,1:4]

#                 PC1          PC2           PC3          PC4
# 10000-28 0.01578677  0.006109059  0.0086136447 -0.010829816
# 10001-7  0.01712065 -0.020603307 -0.0022774886  0.010428056
# 10003-15 0.01163660  0.014874150  0.0088900875  0.012455513
# 10006-25 0.01947500 -0.008470357 -0.0055983740 -0.003819188
# 10008-43 0.01578426  0.003989797 -0.0008580165 -0.017475235

# the matrix x has the principal component score vectors
dim(res.pca$x)
# [1] 1065 1065

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)
# [1] 0.469395583 0.066381867 0.026601245 0.018127898 0.014454832 0.009869478

# Plot cumulative variance explained
plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/eig_values_prcomp.csv")
  
# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# The function princomp() uses the spectral decomposition approach. 
# The functions prcomp() and PCA()[FactoMineR] use the singular value decomposition (SVD).

# Spectral decomposition which examines the covariances / correlations between variables
# Singular value decomposition which examines the covariances / correlations between individuals

###################################################

# Try to create distance matrix with dist()
distance <- dist(prot, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)

# Try mds with k set to max n-1 componenets
out <- cmdscale(distance, k = 1064, eig = TRUE)

GOF <- cmdscale(distance, k = 1064, eig = TRUE)$GOF

coords <- out$points

NewDists <- dist(coords, diag=TRUE, upper=TRUE)

r <- cor(c(distance), c(NewDists))




########################################################################

screen

R

library(tidyverse)
library(ggplot2)
library(readxl)
library(psych)
library(ggcorrplot)
library(cowplot)
library(tidyverse)

# Read in the protein file 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv", check.names = F)

# Isolate just the protein columns of interest for PCA (100 test proteins)
prot <- prot[c(33:132)]

# Try to create distance matrix with dist()
distance <- dist(prot, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)

# Try mds with k set to max n-1 componenets
mds1 <- cmdscale(distance, k = 99)











########################################################################


# Isolate just the protein columns of interest for PCA
prot <- prot[c(33:4267)]

pca <- principal(scale(prot), rotate="none", nfactors=1065)

scores_pca <- pca$scores
variance_pca <- pca$Vaccounted

scores <- as.data.frame(scores_pca)
var1 <- as.data.frame(variance_pca)

var <- var1[1,] # get variance 
cum <- var1[5,] # get cumulative variance

var <- gather(var) # gather
var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
var$num <- as.integer(1:1065)
var$mes <- as.numeric(var$value)

cum <- gather(cum) # gather
cum$num <- as.integer(1:1065)
cum$mes <- as.numeric(cum$value)


var_subset <- var[c(1:500),]

a <- ggplot(var_subset, aes(num, mes, col = col)) +
geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/test_PCA_proteins_eingenvalues_1065_components.pdf", height =20 , width = 20)
a
dev.off()


b <- ggplot(cum, aes(num, mes)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/test_PCA_proteins_cumulative_variance_1065_components.pdf", height =20 , width = 20)
b
dev.off()



# Try prcomp method 
pr <- prcomp(prot, scale = TRUE)







pca <- principal(scale(prot), rotate="none", nfactors=1000)

scores_pca <- pca$scores
variance_pca <- pca$Vaccounted

scores <- as.data.frame(scores_pca)
var1 <- as.data.frame(variance_pca)

var <- var1[1,] # get variance 
cum <- var1[5,] # get cumulative variance

var <- gather(var) # gather
var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
var$num <- as.integer(1:ncol(prot))
var$mes <- as.numeric(var$value)

cum <- gather(cum) # gather
cum$num <- as.integer(1:ncol(prot))
cum$mes <- as.numeric(cum$value)




# Scale and perform PCA
pca <- principal(scale(prot), rotate="none", nfactors=ncol(prot))
# When trying this on 2000 proteins:
# The determinant of the smoothed correlation was zero.
# This means the objective function is not defined for the null model either.
# The Chi square is thus based upon observed correlations.
# In factor.stats, the correlation matrix is singular, an approximation is used

# Warning messages:
# 1: In cor.smooth(model) :
#   Matrix was not positive definite, smoothing was done
# 2: In cor.smooth(r) : Matrix was not positive definite, smoothing was done
# 3: In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   the model inverse times the r matrix is singular, replaced with Identity matrix which means fits are wrong
# 4: In fa.stats(r = r, f = f, phi = phi, n.obs = n.obs, np.obs = np.obs,  :
#   In factor.stats, the correlation matrix is singular, and we could not calculate the beta weights for factor score estimates
# 5: In principal(scale(prot), rotate = "none", nfactors = ncol(prot)) :
#   The matrix is not positive semi-definite, scores found from Structure loadings

scores_pca <- pca$scores
variance_pca <- pca$Vaccounted

scores <- as.data.frame(scores_pca)
var1 <- as.data.frame(variance_pca)

var <- var1[1,] # get variance 
cum <- var1[5,] # get cumulative variance

var <- gather(var) # gather
var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
var$num <- as.integer(1:ncol(prot))
var$mes <- as.numeric(var$value)

cum <- gather(cum) # gather
cum$num <- as.integer(1:ncol(prot))
cum$mes <- as.numeric(cum$value)

# ## Read in seq-id conversion file 
# anno <- read_excel("Y:/Protein_DNAm_Proxies/Annotations_for_reference.xlsx")
# anno <- as.data.frame(anno)
# anno <- anno[,c(1,4,18)]

# ## subset seq-ids
# anno1 = anno[which(anno$SeqId %in% colnames(cor)),] 
# ## match up their order 
# ids = colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] 
# anno1 = anno1[match(ids, anno1$SeqId),]

# ## check they match
# table(anno1$SeqId == colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] )

# ## replace seq-ids with gene names 
# names(anno1)[3] <- "Name"
# colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)
# row.names(cor)[which(row.names(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)


var_subset <- var[c(1:500),]

a <- ggplot(var_subset, aes(num, mes, col = col)) +
geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/test_PCA_proteins_eingenvalues.pdf", height =20 , width = 20)
a
dev.off()


b <- ggplot(cum, aes(num, mes)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/test_PCA_proteins_cumulative_variance.pdf", height =20 , width = 20)
b
dev.off()


# # c <- if(ncol(joint) > 12){  
# c <- ggcorrplot(cor, 
#            hc.order = TRUE)

# c <- ggcorrplot(cor, 
#            hc.order = TRUE,
#            type = "lower")


# c <- c + theme(
#   panel.background = element_rect(fill = "white", colour = "white",
#                                 size = 2, linetype = "solid"),
#   panel.grid.major = element_line(size = 0.5, linetype = 'solid',
#                                 colour = "white"), 
#   panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
#                                 colour = "white")
#   )


# p1 = plot_grid(c,a,b, nrow = 1, labels = c("a", "b", "c"), rel_widths = c(0.85,0.5,0.5))


# pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/test_PCA_2000_proteins.pdf", height =6 , width = 17)
# p1
# dev.off()
