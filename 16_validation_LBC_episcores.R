Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# Validation in two test sets ##############################################
############################################################################################
############################################################################################

# Collation of the LBC1936 projections into STRADL and LBC1921 test sets 
# KORA assessment is performed separately in the KORA projection > STRADL script

library(tidyverse)

# Read in initial holdout correlation results tables with those that passed thresholds (18, 18)
n <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/neuro_pass_new.csv")
i <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/inflam_pass_new.csv")

proxies <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/LBC1921_162_projections.csv")

# Check to make sure all proxies had features in 12 cv first 
length(which(n$Protein %in% colnames(proxies)))
length(which(i$Protein %in% colnames(proxies)))
not <- which(!i$Protein %in% colnames(proxies))
not <- i[not,] # MCP.2 didnt have features in the 12 cv elnets 

############################################################################################

### Testing LBC1936 trained episcores againat LBC1921

############################################################################################

## File of projections to read in for comparisons 
LBC1921 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/LBC1921_162_projections.csv")

# Work out how many had features in 12 cv 
overlap <- which(colnames(LBC1921) %in% n$Protein)
length(overlap) # 19
overlap <- which(colnames(LBC1921) %in% i$Protein)
length(overlap) # 21

### DO CORRELATIONS 
# Read in proteins file (which has already been normalised, residualised and scaled by riccardo)
proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/proteins_162_joint_with_target.csv")
proteins$match <- proteins$Basename
mean <- mean(proteins$age, na.rm = T)
sd <- sd(proteins$age, na.rm = T)

# Read in proxies file 
proxies <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/LBC1921_162_projections.csv")
proxies$match <- proxies$ID

# Subset proxies to only those which have measurements
overlap <- which(colnames(proxies) %in% colnames(proteins))
proxies <- proxies[,overlap]

# Subset proteins to only those which have proxies 
overlap <- which(colnames(proteins) %in% colnames(proxies))
proteins <- proteins[,overlap]

# Match the order of the match identified column between both 
proxies <- proxies[match(proteins$match, proxies$match),]

# Make sure they are identical 
identical(proteins$match, proxies$match)


### FILTER TO UPDATED THRESHOLD LIST

# Load in the updated thresholded list 
neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/neuro_pass_new.csv")
overlap1 <- which(colnames(proteins) %in% neuro$Protein) 
proteins <- proteins[,overlap1] # 19
overlap2 <- which(colnames(proxies) %in% neuro$Protein) 
proxies <- proxies[,overlap2] # 19
# Make results table 
results <- data.frame(Protein = colnames(proxies)[1:19], r = 1:19, p = 1:19, LC = 1:19, UC = 1:19)
# Correlate for each protein against the proxy 
for (i in 1:19){
name <- as.character(colnames(proxies)[i])
cor1 <- cor.test(proteins[,name], proxies[,name])
int <- cor1$conf.int[1:2]
p <- cor1$p.value 
r <- cor1$estimate 
results[i,"r"] <- r
results[i,"p"] <- p
results[i,4:5] <- int
}
results <- results[rev(order(results$r)),]
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_validation_in_LBC1921_neuro.csv", row.names = F)

### PLOT CORRELATIONS AS SCATTER PLOTS FOR PROXIES VS PROTEINS 
# Plot a correlation between projected protein proxies and the normalised protein measurements
library(ggplot2)
library(ggpubr)
library(tidyverse)
plot_list <- list()
for (i in 1:19){
name <- as.character(colnames(proxies)[i])
proxy <- proxies[,name] %>% as.data.frame()
protein <- proteins[,name] %>% as.data.frame()
corrdata <- cbind(proxy, protein)
names(corrdata) <- c("Proxy", "Protein")
p <- ggplot(corrdata, aes(x=Protein, y=Proxy)) +
geom_point(colour = "grey", size = 2) +
geom_smooth(method='lm', colour = "grey") +
theme(axis.title.x=element_text(), axis.text.x=element_text(size=12),     
      axis.title.y=element_text(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = ) + ggtitle(name) +
xlab("Protein measurement") + 
ylab("DNAm projection")
plot_list[[i]] <- p
}
pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_validation_in_LBC1921_neuro_corr_plots.pdf"))
for (i in 1:19) {
    print(plot_list[[i]])
}
dev.off()

############################################################################################

### Testing LBC1936 trained episcores aginst STRADL 

############################################################################################

## File of projections to read in for comparisons 
STRADL <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/STRADL_combined_778_projections.csv")

##File of the residualised proteins 
# resid <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_778_residualised.rds")
resid <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778.csv")

# Make new file match the old file commented out for the proteins 
resid <- resid[-1]
names(resid)[1] <- "Stradl_id"

## Matched olink to soma IDs
total1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/Olink_linked_to_Somascan_IDs.csv")

## Read in STRADL proxy information
proxy <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/STRADL_combined_778_projections.csv", check.names = F)
names(proxy)

## Subset proxy file to the relevant predictors - 41 changed to 47 unique proteins + ID 
proxy_subset <- proxy[,c(1, which(names(proxy) %in% unique(total1$Olink_Protein)))]

## Link up DNAm ids with proteomic IDs 
linker_DNAm = read.table("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_DNAm_target_REM_17April2020.txt", header = T) 

## Subset proteins file to those with DNAm data
proxy_subset = proxy_subset[which(proxy_subset$ID %in% linker_DNAm$DNAm_id),]

proxy_lost = proxy_subset[-which(proxy_subset$ID %in% linker_DNAm$DNAm_id),]

which(linker_DNAm$DNAm_id %in% proxy_lost$ID)

## Subset DNAm data to final protein data 
linker_DNAm = linker_DNAm[which(linker_DNAm$DNAm_id %in% proxy_subset$ID),]

## Make sure both files match up 
id = proxy_subset$ID 
linker_DNAm = linker_DNAm[match(id, linker_DNAm$DNAm_id), ]

## check if they line-up
table(as.character(linker_DNAm$DNAm_id) == as.character(proxy_subset$ID))

## Add in STRADL Sample IDs to proxy_subset file 
proxy_subset$Stradl_id <- linker_DNAm$Stradl_id
proxy_subset <- proxy_subset[,c(48,1:47)]
names(proxy_subset)[2] <- "DNAm_Id"

##### Now, the relevant IDs are the Stradl_ID column in proxy_subset and ID in 
## the prot file 

## Match up IDs in these two files 
prot <- resid
id = proxy_subset$Stradl_id
names(prot)[1] <- "ID"
prot = prot[match(id, prot$ID),]
## check if they line-up
table(as.character(proxy_subset$Stradl_id) == as.character(prot$ID))
#778 - total used 
#810 - total meth generated 

### Metrics - correlations
## Set up column to record correlation coefficient 
total1$r <- 0
total1$p <- 0
total1$LC <- 0
total1$UC <- 0

# Filter to just proteins available in results 
overlap <- which(total1$Olink_Protein %in% colnames(proxy_subset))
total2 <- total1[overlap,]

## Loop through pairs of Olink gene/protein and corresponding Somamer (Uniprot ID) in Somascan data
for(i in 1:nrow(total2)) { 
protein = as.character(total2[i, "Olink_Protein"])
column = as.numeric(total2[i, "Somascan_Column"])
cor1 <- cor.test(proxy_subset[,protein], prot[,column], use = "complete.obs")
int <- cor1$conf.int[1:2]
p <- cor1$p.value 
r <- cor1$estimate 
int <- cor1$conf.int[1:2]
total2[i,"r"] <- r
total2[i,"p"] <- p
total2[i,7:8] <- int
#total2[i,7:8] <- int
} 

### FILTER TO THE UPDATED THRESHOLD LIST OF PROTEINS 
# Read in initial holdout correlation results tables with those that passed thresholds (18, 18)
neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/neuro_pass_new.csv")
inflam <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/inflam_pass_new.csv")
overlap1 <- which(total2$Olink_Protein %in% inflam$Protein) 
inflam_filtered <- total2[overlap1,]
overlap2 <- which(total2$Olink_Protein %in% neuro$Protein)
neuro_filtered <- total2[overlap2,]
inflam_filtered1 <- inflam_filtered[order(inflam_filtered$r, decreasing = TRUE),] 
neuro_filtered <- neuro_filtered[rev(order(neuro_filtered$r)),]
write.csv(inflam_filtered1, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/STRADL_inflam_correlations.csv", row.names = F)
write.csv(neuro_filtered, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/STRADL_neuro_correlations.csv", row.names = F)

### COMBINE WITH EXISTING HOLDOUT RESULTS 
# Read in initial holdout correlation results tables with those that passed thresholds (18, 18)
n <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/neuro_pass_new.csv")
i <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/inflam_pass_new.csv")
# Read in STRADL correlations 
si <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/STRADL_inflam_correlations.csv")
sn <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/STRADL_neuro_correlations.csv")
names(si)[3] <- "Protein"
names(sn)[3] <- "Protein"
library(tidyverse)
inflam <- left_join(i, si, by = "Protein")
neuro <- left_join(n, sn,, by = "Protein")

# Plot a correlation between projected protein proxies and the normalised protein measurements
prot1 <- prot
proxy <- proxy_subset
library(ggplot2)
library(ggpubr)
library(tidyverse)
plot_list <- list()
for(i in 1:47){ 
    protein = as.character(total2[i, "Olink_Protein"])
    column = as.numeric(total2[i, "Somascan_Column"])
    proxy_data <- proxy[,protein] %>% as.data.frame()
    prot_data <- prot1[,column] %>% as.data.frame()
    corrdata <- cbind(proxy_data, prot_data)
    names(corrdata) <- c("Proxy", "Protein")
    p <- ggplot(corrdata, aes(x=Protein, y=Proxy)) +
    geom_point(colour = "grey", size = 2) +
    geom_smooth(method='lm', colour = "grey") +
    theme(axis.title.x=element_text(), axis.text.x=element_text(size=12),     
          axis.title.y=element_text(), axis.text.y=element_text(size=12)) + 
    stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = ) + ggtitle(protein) +
    xlab("Protein measurement") + 
    ylab("DNAm projection")
    plot_list[[i]] <- p
} 

pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_validation_in_STRADL_corr_plots_in_778.pdf"))
for (i in 1:47) {
    print(plot_list[[i]])
}
dev.off()
