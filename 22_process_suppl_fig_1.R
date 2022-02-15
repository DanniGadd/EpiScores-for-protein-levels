Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## SUPPL FIG 1 EPISCORES ###############################
####################################################################################
#################################################################################### 

### CORRELATION PLOTS FOR COVARIATES AND EPISCORES

library(pheatmap)
library(readxl)

# Read in 109 protein scores
d1 <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/files_for_supp_figs/110_protein_scores.csv")
## Remove leading 'X' for numeric seq-ids 
colnames(d1)[28:ncol(d1)] <- gsub("X", "", colnames(d1)[28:ncol(d1)])

## Read in seq-id conversion file
anno <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/somascan_annotation.xlsx")
anno <- as.data.frame(anno)
anno <- anno[,c(1,3)]
anno$Predictor <- gsub("-",".", anno$`Identifier (SOMAScan SeqId or Olink name)`)

## subset seq-ids just to those in 110 proteins
anno1 = anno[which(anno$Predictor %in% colnames(d1)[28:ncol(d1)]),] 
## match up their order 
ids = colnames(d1)[28:ncol(d1)]
anno1 = anno1[match(ids, anno1$Predictor),]
## check they match
table(anno1$Predictor == colnames(d1)[28:ncol(d1)])
## replace seq-ids with gene names 
colnames(d1)[28:ncol(d1)] <- anno1$`Entrez Gene Name`

## Prepare correlation matrix 
row.names(d1) <- d1$Sample_Name
d1$Sample_Name <- NULL 
d1 = d1[,-grep("IL.12", names(d1))]

## Make correlation matrix
cor <- cor(d1)

## Make annotation df - olink + somalogic proteins - matching order as they appear in corr matrix
annotation <- data.frame(Platform = factor(1:109, labels = c(rep("Olink",25), rep("Somalogic", 84))))

## make annotation colors list and dataframe
Platform = c("#CBCE91FF", "#76528BFF")
names(Platform) = c("Olink", "Somalogic")
## have to make unique as some proteins/gene names are present in olink and somalogic e.g. MMP1
colnames(cor) <- make.unique(colnames(cor))
row.names(cor) <- make.unique(row.names(cor))

colnames(cor)[grep("\\.1", colnames(cor))]

colnames(cor)[colnames(cor) %in% "NTRK3.1"] <- "NTRK3.2"
colnames(cor)[colnames(cor) %in% "NTRK3"] <- "NTRK3.1"

colnames(cor)[colnames(cor) %in% "GZMA.1"] <- "GZMA.2"
colnames(cor)[colnames(cor) %in% "GZMA"] <- "GZMA.1"

colnames(cor)[colnames(cor) %in% "CXCL11.1"] <- "CXCL11.2"
colnames(cor)[colnames(cor) %in% "CXCL11"] <- "CXCL11.1"

colnames(cor)[colnames(cor) %in% "CXCL10.1"] <- "CXCL10.2"
colnames(cor)[colnames(cor) %in% "CXCL10"] <- "CXCL10.1"

colnames(cor)[colnames(cor) %in% "CLEC11A.1"] <- "CLEC11A.2"
colnames(cor)[colnames(cor) %in% "CLEC11A"] <- "CLEC11A.1"

row.names(cor)[row.names(cor) %in% "NTRK3.1"] <- "NTRK3.2"
row.names(cor)[row.names(cor) %in% "NTRK3"] <- "NTRK3.1"

row.names(cor)[row.names(cor) %in% "GZMA.1"] <- "GZMA.2"
row.names(cor)[row.names(cor) %in% "GZMA"] <- "GZMA.1"

row.names(cor)[row.names(cor) %in% "CXCL11.1"] <- "CXCL11.2"
row.names(cor)[row.names(cor) %in% "CXCL11"] <- "CXCL11.1"

row.names(cor)[row.names(cor) %in% "CXCL10.1"] <- "CXCL10.2"
row.names(cor)[row.names(cor) %in% "CXCL10"] <- "CXCL10.1"

row.names(cor)[row.names(cor) %in% "CLEC11A.1"] <- "CLEC11A.2"
row.names(cor)[row.names(cor) %in% "CLEC11A"] <- "CLEC11A.1"

row.names(annotation) <- colnames(cor)
ann_colors = list(Platform = Platform)


## Store resulting pheatmap in pdf 
pdf("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Supp_Figs/Supp_Fig1.pdf", height = 9.4, width = 9.4)
pheatmap(cor,  border_color = NA, treeheight_row= 0, treeheight_col =  0, fontsize =6.3, legend = T, annotation_col = annotation, annotation_colors = ann_colors)
dev.off()

library(readxl)

# Read in 110 protein scores
d1 <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/files_for_supp_figs/110_protein_scores.csv")
## Remove leading 'X' for numeric seq-ids 
colnames(d1)[28:ncol(d1)] <- gsub("X", "", colnames(d1)[28:ncol(d1)])
