Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

######################################################################################
######################################################################################
###################### PCA ON VALIDATED EPISCORES ######################################
######################################################################################
######################################################################################

# PCA analysis to understand whether we have multiple episcore signals for those passing
# the validation step.

library(tidyverse)

### SOMA LOAD IN AND SELECT EPISCORES FOR COX MODELS
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/GS_combined_9537_projections.csv")
names(prot)[1] <- "XID"
colnames(prot) <- sub('.', '', colnames(prot))
#prot = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/olink_files/GS_combined_9537_projections.csv")
# read in proteins taken forward from validation test steps 
results <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Correlations_020221.csv")
results <- results %>% filter(r > 0.1)
results <- results %>% filter(p < 0.05)
list3 <- results$SeqId %>% as.data.frame()
names(list3)[1] <- "Protein"
list3$Protein <- str_replace_all(list3$Protein, "-", ".")
#list3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/olink_files/List_of_proteins_for_cox_models_220121.csv")
ID <- prot$ID
over <- which(colnames(prot) %in% list3$Protein)
pre <- prot[over]
pred <- cbind(ID, pre)
soma <- pred

### OLINK LOAD IN AND SELECT EPISCORES FOR COX MODELS ()
prot = read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_projections/GS_combined_9537_projections.csv")
list3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Choosing_list_of_top_proxies/List_of_proteins_for_cox_models_220121.csv")
ID <- prot$ID
over <- which(colnames(prot) %in% list3$Protein)
pre <- prot[over]
pred <- cbind(ID, pre)
olink <- pred

### JOIN UP EPISCORES 
# first match the olink to the soma IDs 
people <- soma$ID
matched <- olink[match(people, olink$ID),]
join <- cbind(soma, matched)


### PCA on all the proteins going into cox models 
# Check no NA values 
print(table(is.na(joint))) # there are none
# Load PSYCH 
#install.packages("psych")
library(psych)
joint <- join[-1] # remove the ID column
joint <- joint[-136] # remove the second ID column
# > dim(joint)
# [1] 9537  165

# Scale data and get PC scores
scores_pca <- principal(scale(joint), rotate="none", nfactors=40)$scores
variance_pca <- principal(scale(joint), rotate="none", nfactors=40)$Vaccounted
pca <- principal(scale(joint), rotate="none", nfactors=40)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/PCA_on_val/PCA.pdf", width = 10, height = 8)
par(mfrow = c(1,2))
plot(pca$Vaccounted[1,], pch = 16, lty = 2)
plot(pca$Vaccounted[5,], pch = 16, lty = 2)
dev.off()
scores <- as.data.frame(scores_pca)
var <- as.data.frame(variance_pca)
write.csv(scores, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/PCA_on_val/scores_PCA.csv")
write.csv(var, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/PCA_on_val/var_PCA.csv")
scores1 <- as.data.frame(scores_pca)
ggplot(scores, aes(x=PC1, y=PC2)) +
	geom_point()
d <- principal(scale(joint), rotate="none", nfactors=40)
data <- cbind(joint, scores[c(1:2)])
library(ggplot2)
p <- ggplot(data, aes(PC1, PC2, col = "Protein")) +
	stat_ellipse(geom = "polygon", col = "black", alpha = 0.5)+
	geom_point(shape = 21, col = "black")















































