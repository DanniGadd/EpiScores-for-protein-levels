Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# KORA projections #########################################################
############################################################################################
############################################################################################

library(readxl)
library(tidyverse)

# Read in the stradl proteins file which has now been matched to the annotation file with SeqIDs
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/260121_Matching_cohorts/STRADL_778_residualised_matched.csv", check.names = F)

# Weight file 
cpgs <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Weights_230221/weights_KORA.csv")
un <- unique(cpgs$Predictor) %>% as.data.frame() # 793 unique predictors 
names(un) <- "Protein"

# Create matchable ID 
cpgs$ID<- gsub("\\:.*","", cpgs$Predictor)
cpgs$ID <- gsub("\\_.*","", cpgs$ID)
cpgs <- cpgs[c(1,2,4)]
names(cpgs)[3] <- "Predictor"

# Remove intercept
intercepts <- which(cpgs$CpG_site %in% "(Intercept)")
cpgs <- cpgs[-intercepts,]

# This set should be the cpgs now ready to run predictor scores 
# I will first check against the stradl seqIDs and subset based on those that are available across cohorts
overlap <- which(cpgs$Predictor %in% colnames(prot))
cpgs <- cpgs[overlap,]
predictors <- cpgs$Predictor %>% unique()
length(predictors) # 793

# Create a version of the cpg weights which has the uniprot ID linked and the name of the protein 
# Annotation file merge into results 
anno <- read_excel("/Cluster_Filespace/Marioni_Group/STRADL_proteins/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- anno %>% as.data.frame()
test <- anno[c(1:4,13,18,20)]
cpgs2 <- cpgs
names(cpgs2)[3] <- "SeqId"
test2 <- left_join(cpgs2, test, by = "SeqId")

# remove the coefs which = 0 
ov <- which(test2$Coefficient == "0")
test3 <- test2[-ov,]
test3$Panel <- "SomaScan"
test3$Training_cohort <- "KORA"

# count unique proteins with features 
un <- test3$SeqId %>% unique() # 480 
# fiter to those in final set (84)
list <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/list_of_cox_proteins.csv")
ov <- which(test3$SeqId %in% list)
test4 <- test3[ov,]
length(test4$SeqId %>% unique()) # 84
test4 <- test4[c(1:4,6,10,11)]
write.csv(test4, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/supplementary_KORA_cpgs_for_predictors_all_020321.csv", row.names = F)

### ADD protein gene naming to suppl table with predictor weights 
library(tidyverse)
library(readxl)
file <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results/Linker_predictor_weights_file.xlsx")
file <- file %>% as.data.frame()
names(file)[3] <- "SeqId"
linker <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/somascan_annotation.xlsx")
linker <- linker %>% as.data.frame()
names(linker)[1] <- "SeqId"
file2 <- left_join(file, linker, by = "SeqId")
# file2 <- file2[c(1,2,6,7,9,8,3,5,4)]
write.csv(file2, "Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results/Linker_file_protein_naming_corrected.csv", row.names = F)
file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/Linker_file_protein_naming_corrected.csv")
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
names(anno)[4] <- "CpG.site"
anno <- anno[c("CpG.site", "chr", "pos", "UCSC_RefGene_Name")]
join <- left_join(file, anno, by = "CpG.site")
join$UCSC_RefGene_Name <- str_split(join$UCSC_RefGene_Name, ";")
join$UCSC_RefGene_Name <- lapply(join$UCSC_RefGene_Name, unique)
join$UCSC_RefGene_Name <- vapply(join$UCSC_RefGene_Name, paste, collapse = ", ", character(1L))
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/Linker_file_protein_naming_corrected_with_annotation.csv", row.names = F)


############################################################################################

### Project KORA weights into STRADL 

############################################################################################

# Read in DNAm data for STRADL in 778
meth2 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/W2_prep_479.rds")
meth3 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/W3_prep_299.rds")

# Project weights into STRADL wave 2
coef <- t(meth2)

loop = unique(cpgs$Predictor)
out <- data.frame()

#i <- "loop[1]"
i <- "5476-66"
5400-52
5476-66 

for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 

out_w2 <- out 

# Project weights into STRADL wave 3
coef <- t(meth3)

loop = unique(cpgs$Predictor)
out <- data.frame()

for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 

out_w3 <- out 

out <- rbind(out_w2, out_w3)

dim(out)
# [1] 778 793

# Save a copy of the projected proteins in STRADL 
write.csv(out, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/KORA_projections_into_STRADL_778.csv")


## Correlate the stradl proteins from our own dataset against KORA scores
# We will need to join in the DNAm ID from the target file 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_DNAm_target_REM_17April2020.txt")
IDs <- target[c(1,3)]
library(tidyverse)
prot2 <- left_join(prot, IDs, by = "Stradl_id")
prot2 <- prot2[match(row.names(out), prot2$DNAm_id),]
#colnames(prot2) <- gsub("\\..*","", colnames(prot2))
identical(row.names(out), prot2$DNAm_id) # TRUE

# Correlate proxy
dim <- dim(out)[2]
result <- data.frame(SeqId = 1:dim, r = 1:dim, p = 1:dim, LC = 1:dim, UC = 1:dim)
for(i in 1:dim){
  seqID <- colnames(out)[i]
  cor <- cor.test(prot2[,seqID], out[,seqID])
  int <- cor$conf.int[1:2]
  result[i,1] <- seqID
  result[i,2] <- cor$estimate
  result[i,3] <- cor$p.value
  result[i,4:5] <- int
}
result <- result[order(result$r, decreasing = TRUE),] 


# Create a version with NAs removed to see how many didnt generate scores 
res <- na.omit(result) 
dim(res) # 480 generated scores 

# Annotation file merge into results 
anno <- read_excel("/Cluster_Filespace/Marioni_Group/STRADL_proteins/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
test <- left_join(res, anno, by = "SeqId")
test <- test[c(1:8,17,22,24)]
write.csv(test, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/KORA_correlations_in_STRADL_all.csv", row.names = F)

# Now look at which pass our threshold 
pass <- test %>% filter(r > 0.1) # 84
pass <- pass %>% filter(p < 0.05) # 84
write.csv(pass, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/KORA_correlations_in_STRADL_filtered.csv", row.names = F)

### Read file back in to test FDR 
library(readxl)
library(tidyverse)
file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/KORA_correlations_in_STRADL_all.csv")
file$adj <- p.adjust(file$p, method = "BH")
filt <- file %>% filter(adj < 0.05)

##############################################################################

### PLOT CORRELATIONS AS SCATTER PLOTS FOR EPISCORES VS PROTEINS 

# Plot a correlation between projected protein proxies and the normalised protein measurements
proxies <- out
ov <- which(colnames(proxies) %in% res$SeqId) # filter to do plots for the 480 with scores 
proxies <- proxies[,ov]

proteins <- prot2

library(ggplot2)
library(ggpubr)
library(tidyverse)

plot_list <- list()


for (i in 1:480){
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


pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/KORA_validation_in_STRADL_corr_plots_in_778.pdf"))
for (i in 1:480) {
    print(plot_list[[i]])
}
dev.off()


# Repeat but do for the 84 which were robust
proxies <- out
ov <- which(colnames(proxies) %in% pass$SeqId) # filter to do plots for the 84 
proxies <- proxies[,ov]

proteins <- prot2

library(ggplot2)
library(ggpubr)
library(tidyverse)

plot_list <- list()


for (i in 1:84){
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


pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/KORA_validation_in_STRADL_corr_plots_in_778_84_significant.pdf"))
for (i in 1:84) {
    print(plot_list[[i]])
}
dev.off()








