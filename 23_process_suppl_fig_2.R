Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## SUPPL FIG 2 #########################################
####################################################################################
#################################################################################### 

### CORRELATION PLOTS FOR COVARIATES AND EPISCORES

library(pheatmap)
library(readxl)

# Read in 78 proxies with associations with incident disease + covariates

d1 <- read.csv("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/files_for_supp_figs/78_protein_scores.csv")

grimage_w1 <- read.csv("U:/Generation_Scotland_data/DNAmAge/STRADL_w1_DNAmAge_combined_agemonths.csv")
grimage_w3 <- read.csv("U:/Generation_Scotland_data/DNAmAge/STRADL_w3_DNAmAge_combined.csv")

grimage_w1$SampleID <- gsub("X", "", grimage_w1$SampleID)
grimage_w3$SampleID <- gsub("X", "", grimage_w3$SampleID)

id_w1 <- read.csv("U:/Generation_Scotland_data/DNAmAge/stradl-samples-5101.csv")
id_w3 <- read.csv("U:/Generation_Scotland_data/DNAmAge/samplesheet.final.csv")

id_w1$SampleID <- paste(id_w1$Sentrix_ID, id_w1$Sentrix_Position, sep = "_")
id_w3$SampleID <- paste(id_w3$Sentrix_ID, id_w3$Sentrix_Position, sep = "_")

id_w1 <- id_w1[,c("SampleID", "Sample_Name")]
id_w3 <- id_w3[,c("SampleID", "Sample_Name")]

grimage_w1 <- merge(grimage_w1, id_w1, by = "SampleID")
grimage_w3 <- merge(grimage_w3, id_w3, by = "SampleID")

grimage_w1 <- grimage_w1[,c("Sample_Name", "AgeAccelGrim")]
grimage_w3 <- grimage_w3[,c("Sample_Name", "AgeAccelGrim")]

grimage <- rbind(grimage_w1, grimage_w3)

d1 <- merge(grimage, d1, by ="Sample_Name")

## remove leading 'X' from numeric seq-ids
colnames(d1)[34:ncol(d1)] <- gsub("X", "", colnames(d1)[34:ncol(d1)])


## Read in seq-id conversion file 
anno <- read.csv("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
anno <- anno[,c(1,4)]
anno$SeqId <- gsub("-",".", anno$Predictor)

## subset seq-ids just to those in 78 proteins
anno1 = anno[which(anno$SeqId %in% colnames(d1)[34:ncol(d1)]),] 
## match up their order 
ids = colnames(d1)[34:ncol(d1)]
anno1 = anno1[match(ids, anno1$SeqId),]
## check they match
table(anno1$SeqId == colnames(d1)[34:ncol(d1)])
## replace seq-ids with gene names 
colnames(d1)[34:ncol(d1)] <- as.character(anno1$Name)


## Extract covariates of interest 
life <- d1[,c("Age", "units", "smokingScore", "simd", "EA", "bmi", "AgeAccelGrim", "CD8T", "CD4T", "Bcell", "NK", "Mono", "Gran")]
names(life) <- c("Age", "Alcohol", "EpiSmoker", "SIMD", "EA", "BMI", "GrimAge", "CD8T", "CD4T", "Bcell", "NK", "Mono", "Gran")
lifename <- names(life)

## Extract protxies 
list <- d1[,c(15:ncol(d1))]

## Set up matrix to loop through identfying correlations between proxies + covariates
mat1 <- matrix(nrow = length(list), ncol = length(lifename)) 

## Loop for correlations 
for(i in 1:ncol(list)){
  for(j in 1:ncol(life)){ 
  mat1[i,j] <- cor(list[,i], life[,j], use = "complete.obs")
  }
} 

mat1 = as.data.frame(mat1)

## Input names for correlations 
colnames(mat1) <- lifename
row.names(mat1) <- make.unique(names(d1)[15:ncol(d1)])

row.names(mat1)[row.names(mat1) %in% "NTRK3"] <- "NTRK3.1"


## Get annotations ready for pheatmap 
annotationrow <- data.frame(Platform = factor(1:78, labels = c(rep("Olink",19), rep("Somalogic", 59))))
annotationcol <- data.frame(Type = factor(1:13, labels = c(rep("Life",7), rep("White Blood Cell", 6))))


## make annotation colors list and dataframe

Platform = c("#CBCE91FF", "#76528BFF")
names(Platform) <- c("Olink", "Somalogic")
Type = c("#3eb489", "#FFB6C1")
names(Type) <- c("Life", "White Blood Cell")
ann_colors = list(Platform=Platform, Type = Type)


## have to make unique as some proteins/gene names are present in olink and somalogic e.g. NTRK3
rownames(annotationrow) <- make.unique(colnames(list))
rownames(annotationcol) <- make.unique(colnames(life)) ## wont use this in plot but keeping it in case

## Make nice blue, white, red gradient for strength of correlations 
myColor <- colorRampPalette(c("blue", "white", "red"))(300)

## Store plot in pdf 

pdf("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Supp_Figs/Supp_Fig2.pdf", height = 8.5, width =10.7)
pheatmap(mat1, color = myColor, display_numbers = F, cluster_rows = F, cluster_cols = F, breaks = seq(-1,1,length.out = 301),fontsize_col = 10.5, fontsize_row = 7.5, angle_col = 0, annotation_row = annotationrow, annotation_names_row = F, annotation_colors = ann_colors)
dev.off()




