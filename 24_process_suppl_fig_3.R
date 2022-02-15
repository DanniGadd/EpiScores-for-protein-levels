Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## EPISCORE PLOTS ######################################
####################################################################################
#################################################################################### 

### PCA ON MULTIPLE PROXY PREDICTORS FOR DISEASE

library(readxl)
library(psych)
library(ggcorrplot)
library(cowplot)
library(tidyverse)

# One plot for each disease 

####################################################################################

### GET d1 proteins file 

####################################################################################

d1 <- read.csv("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/files_for_supp_figs/78_protein_scores.csv")
colnames(d1)[33:ncol(d1)] <- gsub("X", "", colnames(d1)[33:ncol(d1)])

## Read in Disease Data
diseases <-read.csv("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/files_for_supp_figs/FullyAdjusted_Disease_Proteins.csv")
diseases$Protein <- gsub("-", ".", diseases$Protein)

loop = list.files("U:/Protein_DNAm_Proxies/Cox_Cases", ".csv")

setwd("U:/Protein_DNAm_Proxies/Cox_Cases/")
loop1 = loop[c(10,12,15,7,14,11,8,6)]

list_diseases <- gsub("_all.*", "", loop1)

list_plots <- list()
for(i in 1:6){ 
  
tmp = read.csv(loop1[i])
dis <- list_diseases[i]

for_pca <- diseases[which(diseases$Disease %in% dis),]

joint <- d1[,c(1,which(colnames(d1) %in% unique(for_pca$Protein)))]

joint <- joint[which(joint$Sample_Name %in% tmp$Sample_Name), ]
joint$Sample_Name <- NULL

# Scale data and get PC scores
scores_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$scores
variance_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$Vaccounted
pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))

scores <- as.data.frame(scores_pca)
var1 <- as.data.frame(variance_pca)

var <- var1[1,] # get variance 
cum <- var1[5,] # get cumulative variance

var <- gather(var) # gather
var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
var$num <- as.integer(1:ncol(joint))
var$mes <- as.numeric(var$value)

cum <- gather(cum) # gather
cum$num <- as.integer(1:ncol(joint))
cum$mes <- as.numeric(cum$value)

cor = cor(joint)

## Read in seq-id conversion file 
anno <- read.csv("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
anno <- anno[,c(1,4)]
anno$SeqId <- gsub("-",".", anno$Predictor)

## subset seq-ids just to those in 78 proteins
anno1 = anno[which(anno$SeqId %in% colnames(cor)),] 
## match up their order 
ids = colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] 
anno1 = anno1[match(ids, anno1$SeqId),]
## check they match
table(anno1$SeqId == colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] )
## replace seq-ids with gene names 
colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)
row.names(cor)[which(row.names(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)

a <- if(ncol(joint) > 12) { ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
} else { 
  
  ggplot(var, aes(num, mes, col = col)) +
    geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
    xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
    scale_color_manual(values=c("#E69F00", "#999999")) + scale_x_continuous(n.breaks = ncol(joint))
  
  
}


b <- if(ncol(joint) > 12) { ggplot(cum, aes(num, mes)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 
} else if(ncol(joint) > 4) { ggplot(cum, aes(num, mes)) +
  geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
  xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) +  scale_x_continuous(n.breaks = ncol(joint))} else { 
    ggplot(cum, aes(num, mes)) +
      geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
      xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 
    
    }

c <- if(ncol(joint) > 12){  ggcorrplot(cor, 
           hc.order = TRUE, 
           lab = TRUE,
           type = "lower",
           lab_size = 1.5,
           colors = c("blue", "white", "red")) +labs(title = list_diseases[i]) + 
  theme(legend.title = element_blank(), text = element_text(size = 11), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
} else { 
  ggcorrplot(cor, 
             hc.order = TRUE, 
             lab = TRUE,
             type = "lower",
             lab_size = 2.5,
             colors = c("blue", "white", "red")) +labs(title = list_diseases[i]) + 
    theme(legend.title = element_blank(), text = element_text(size = 11), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
  
  } 

p1 = plot_grid(c,a,b, nrow = 1, labels = c("a", "b", "c"), rel_widths = c(0.85,0.5,0.5))
list_plots[[i]] <- p1

print(i)
} 


pdf("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Supp_Figs/Fig3_Diseases1-3.pdf", height = 11, width = 10.5)
plot_grid(list_plots[[1]],list_plots[[2]],list_plots[[3]], nrow = 3)
dev.off()

pdf("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Supp_Figs/Fig3_Diseases4-6.pdf", height = 11, width = 10.5)
plot_grid(list_plots[[4]],list_plots[[5]],list_plots[[6]], nrow = 3)
dev.off()






##### DISEASES WITH LARGE NUMBERS OF ASSOCATIONS #############

list_plots2 <- list()
for(i in 7:8){ 
  
  tmp = read.csv(loop1[i])
  dis <- list_diseases[i]
  
  for_pca <- diseases[which(diseases$Disease %in% dis),]
  
  joint <- d1[,c(1,which(colnames(d1) %in% unique(for_pca$Protein)))]
  
  joint <- joint[which(joint$Sample_Name %in% tmp$Sample_Name), ]
  joint$Sample_Name <- NULL
  
  # Scale data and get PC scores
  scores_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$scores
  variance_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$Vaccounted
  pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))
  
  scores <- as.data.frame(scores_pca)
  var1 <- as.data.frame(variance_pca)
  
  var <- var1[1,] # get variance 
  cum <- var1[5,] # get cumulative variance
  
  var <- gather(var) # gather
  var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
  var$num <- as.integer(1:ncol(joint))
  var$mes <- as.numeric(var$value)
  
  cum <- gather(cum) # gather
  cum$num <- as.integer(1:ncol(joint))
  cum$mes <- as.numeric(cum$value)
  
  cor = cor(joint)
  
  ## Read in seq-id conversion file 
  anno <- read.csv("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
  anno <- anno[,c(1,4)]
  anno$SeqId <- gsub("-",".", anno$Predictor)
  
  ## subset seq-ids just to those in 78 proteins
  anno1 = anno[which(anno$SeqId %in% colnames(cor)),] 
  ## match up their order 
  ids = colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] 
  anno1 = anno1[match(ids, anno1$SeqId),]
  ## check they match
  table(anno1$SeqId == colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] )
  ## replace seq-ids with gene names 
  colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)
  row.names(cor)[which(row.names(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)
  
  row.names(cor) <- make.unique(row.names(cor))
  colnames(cor) <- make.unique(colnames(cor))
  
  a <-  ggplot(var, aes(num, mes, col = col)) +
      geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
      xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
      scale_color_manual(values=c("#E69F00", "#999999")) 
 

  
  b <- ggplot(cum, aes(num, mes)) +
      geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
      xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 

      
  
  c <- ggcorrplot(cor, 
                                         hc.order = TRUE, 
                                         lab = F,
                                         type = "lower",
                                         colors = c("blue", "white", "red")) +labs(title = list_diseases[i]) + 
      theme(legend.title = element_blank(), text = element_text(size = 5.5), axis.text.x=element_text(size=5.5), axis.text.y=element_text(size=5.5))

  
  p1 = plot_grid(c,a,b, nrow = 1, labels = c("a", "b", "c"), rel_widths = c(0.85,0.5,0.5), rel_heights = c(0.7,0.5,0.5))
  list_plots2[[i]] <- p1
  
  print(i)
} 


pdf("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Supp_Figs/Fig3_Diseases7-9.pdf", height = 8, width = 10)
plot_grid(list_plots2[[7]],list_plots2[[8]], nrow = 2)
dev.off()













