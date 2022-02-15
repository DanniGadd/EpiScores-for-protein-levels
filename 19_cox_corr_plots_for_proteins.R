Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## CORR PLOTS ##########################################
####################################################################################
####################################################################################

# Correlation plots for top episcores/proxies against key variables 

library(tidyverse)

### SOMA LOAD IN AND SELECT PROXIES FOR COX MODELS (135)
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

### OLINK LOAD IN AND SELECT PROXIES FOR COX MODELS ()
prot = read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_projections/GS_combined_9537_projections.csv")
# read in proteins taken forward from validation test steps 
list3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Choosing_list_of_top_proxies/List_of_proteins_for_cox_models_220121.csv")
ID <- prot$ID
over <- which(colnames(prot) %in% list3$Protein)
pre <- prot[over]
pred <- cbind(ID, pre)
olink <- pred

### JOIN UP PROXIES 

# first match the olink to the soma IDs 
people <- soma$ID
matched <- olink[match(people, olink$ID),]
join <- cbind(soma, matched)
join <- join[-137]

join <- join[-1]


####################################################################################

### Corr plots 

####################################################################################

data <- join 

cor <- cor(data)

library(ggcorrplot)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/corr_plot_for_proteins_164.pdf",width = 40, height = 40)
ggcorrplot(cor, hc.order = TRUE, type = "lower",
     outline.col = "white")
dev.off()


# Restrict to those which had disease associations 

results <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/Joint_results_with_olink/results_cox_040221.csv")

proteins <- results$SeqId %>% unique() # 82

data2 <- join 

data2 <- data2[which(colnames(data2) %in% proteins)]


cor <- cor(data2)

library(ggcorrplot)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/corr_plot_for_proteins_associated_with_disease_82.pdf",width = 30, height = 30)
ggcorrplot(cor, hc.order = TRUE, type = "lower",
     outline.col = "white")
dev.off()




####################################################################################

### Grimage vs proxies 

####################################################################################

### Reload the following and join the join object to d1 to get age accel scores 


library(tidyverse)

### SOMA LOAD IN AND SELECT PROXIES FOR COX MODELS (135)

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

### OLINK LOAD IN AND SELECT PROXIES FOR COX MODELS ()

prot = read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_projections/GS_combined_9537_projections.csv")

# read in proteins taken forward from validation test steps 

list3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Choosing_list_of_top_proxies/List_of_proteins_for_cox_models_220121.csv")

ID <- prot$ID
over <- which(colnames(prot) %in% list3$Protein)
pre <- prot[over]
pred <- cbind(ID, pre)

olink <- pred

### JOIN UP PROXIES 

# first match the olink to the soma IDs 
people <- soma$ID
matched <- olink[match(people, olink$ID),]
join <- cbind(soma, matched)

join <- join[-137]


# d1 file

d1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/d1.csv")

cols <- d1[c(1,2,3,4,25,108,109)]
names(cols)[2] <- "ID"

joint <- left_join(join, cols, by = "ID")

names <- colnames(joint)[2:165]

### Plot proxies vs grimage accel 

library(ggpubr)

# Set dataset 
data <- joint
data$AgeAccelGrim <- as.numeric(data$AgeAccelGrim)

# Plot
plot_list <- list()

for(i in names){
      id <- i
      interest <- data[c("AgeAccelGrim", id)]
      names(interest)[2] <- "protein"
      p <- ggplot(interest, aes(x=AgeAccelGrim, y=protein)) +
        geom_point() + geom_smooth(method = 'lm', color = "grey") + theme_minimal() + ggtitle(i) + ylab(i) + xlab("GrimAge acceleration")
        stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7)
      plot_list[[i]] <- p 
}

pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/grimage_accel_vs_proxies.pdf"))
for (i in 1:164) {
    print(plot_list[[i]])
}
dev.off()


### Record the correlations for the plots

# Make results table 
results <- data.frame(Proxy = names, r = 1:29, p = 1:29, LC = 1:29, UC = 1:29)

# Correlate for each protein against the proxy 

for (i in 1:29){
name <- names[i]
cor1 <- cor.test(data[,name], data$AgeAccelGrim.y)
int <- cor1$conf.int[1:2]
p <- cor1$p.value 
r <- cor1$estimate 
results[i,"r"] <- r
results[i,"p"] <- p
results[i,4:5] <- int
}

# Order by r
results <- results[rev(order(results$r)),]

write.csv(results, "/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Grimage/grimage_accel_vs_proxies_results_table.csv", row.names = F)


####################################################################################

### Age vs proxies

####################################################################################

### Plot proxies vs grimage accel 

library(ggpubr)

names <- read.csv("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Choosing_list_of_top_proxies/List_of_proteins_for_cox_models_220121.csv")
names <- names$Protein %>% as.character()

# Set dataset 
data <- join
data$Age <- as.numeric(data$Age)

# Plot
plot_list <- list()

for(i in names){
      id <- i
      interest <- data[c("Age", id)]
      names(interest)[2] <- "protein"
      p <- ggplot(interest, aes(x=Age, y=protein)) +
        geom_point() + geom_smooth(method = 'lm', color = "grey") + theme_minimal() + ggtitle(i) + ylab(i) + xlab("Age")
        stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7)
      plot_list[[i]] <- p 
}

pdf(file = paste0("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Grimage/age_vs_proxies.pdf"))
for (i in 1:29) {
    print(plot_list[[i]])
}
dev.off()


### Record the correlations for the plots

# Make results table 
results <- data.frame(Proxy = names, r = 1:29, p = 1:29, LC = 1:29, UC = 1:29)

# Correlate for each protein against the proxy 

for (i in 1:29){
name <- names[i]
cor1 <- cor.test(data[,name], data$Age)
int <- cor1$conf.int[1:2]
p <- cor1$p.value 
r <- cor1$estimate 
results[i,"r"] <- r
results[i,"p"] <- p
results[i,4:5] <- int
}

# Order by r
results <- results[rev(order(results$r)),]

write.csv(results, "/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Grimage/age_vs_proxies_results_table.csv", row.names = F)


#################################################################################################################

# Plot proxy levels by Sex in the GS population

#################################################################################################################

names <- read.csv("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Choosing_list_of_top_proxies/List_of_proteins_for_cox_models_220121.csv")
names <- names$Protein %>% as.character()

# Set dataset 
data$Female <- as.character(data$Female)

# Plot 
plot_list <- list()

for(i in names){
      id <- i
      interest <- data[c("Female", id)]
      names(interest)[2] <- "protein"
      p <- ggplot(interest, aes(x=Female, y=protein)) +
        geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
        geom_boxplot(width=0.1) + theme_minimal() + ggtitle(i)

      plot_list[[i]] <- p 
}

pdf(file = paste0("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Grimage/proxies_by_sex.pdf"))
for (i in 1:29) {
    print(plot_list[[i]])
}
dev.off()



#################################################################################################################

# Plot proxy levels by Set in the GS population

#################################################################################################################


# Set dataset 
data$Set <- as.character(data$Set)

# Plot each protein proxy by set 1 and 2 as violin boxplots first 
plot_list <- list()

for(i in names){
      id <- i
      interest <- data[c("Set", id)]
      names(interest)[2] <- "protein"
      p <- ggplot(interest, aes(x=Set, y=protein)) +
        geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
        geom_boxplot(width=0.1) + theme_minimal() + ggtitle(i)

      plot_list[[i]] <- p 
}

pdf(file = paste0("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Grimage/proxies_by_set.pdf"))
for (i in 1:29) {
    print(plot_list[[i]])
}
dev.off()




#################################################################################################################

# Disease specific plots 

#################################################################################################################

# read in the cox tables for each trait

AD <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/AD_extracted_cases_over_60_years_old.csv")
bowel <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Bowel.Cancer_all_cases.csv")
stroke <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Stroke_all_cases.csv")
RA <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/RA_all_cases.csv")
pain <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Pain_all_cases.csv")
lung <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Lung.Cancer_all_cases.csv")
IHD <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/IHD_all_cases.csv")
diab <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Diabetes_all_cases.csv")
depr <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Depression_all_cases.csv")
IBD <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/IBD_all_cases.csv")
COPD <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/COPD_all_cases.csv")
breast <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Cox_tte_tables_for_12_traits/Breast.Cancer_all_cases.csv")

#################################################################################################################

# Disease specific plots - plot by age in each case/control case (f.adj and basic versions for comparison as cases change)

#################################################################################################################

### BASIC MODELS 

data <- d1 

my.list = list(AD, bowel, stroke, RA, pain, lung, IHD, diab, depr, IBD, COPD, breast)
names <- c("AD", "bowel", "stroke", "RA", "pain", "lung", "IHD", "diab", "depr", "IBD", "COPD", "breast")

# Plot each disease by cases and controls vs age 
plot_list <- list()

for(i in 1:12){
  disease <- my.list[[i]]
  join <- left_join(disease, data, by = "Sample_Name")
  # add status identifier 
  join <- join %>% mutate(Status = case_when(
  Event == "1" ~ "Cases",
  Event == "0" ~ "Controls"))
  # remove NA event rows
  overlap <- which(!join$Event == "NA")
  join <- join[overlap,]
  # get relevant info for plot 
  interest <- join[c("Status", "age_at_event")]
  # violin plot 
  p <- ggplot(interest, aes(x=Status, y=age_at_event)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ylab("Age") + ggtitle(names[[i]]) + xlab(names[[i]])
  plot_list[[i]] <- p
}
  
pdf(file = paste0("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Grimage/All_disease_cases_and_controls_by_age_basic_model.pdf"))
for (i in 1:12) {
    print(plot_list[[i]])
}
dev.off()






#################################################################################################################

### PLOTS FOR EACH DISEASE AND PROTEIN (f.adj and basic versions for comparison as cases change)

#################################################################################################################


### FULLY ADJUSTED 

data <- d1 

my.list = list(AD, bowel, stroke, RA, pain, lung, IHD, diab, depr, IBD, COPD, breast)
names <- c("AD", "bowel", "stroke", "RA", "pain", "lung", "IHD", "diab", "depr", "IBD", "COPD", "breast")


names <- read.csv("/Volumes/marioni-lab/Danni/LBC_proteins_Jan2021/Choosing_list_of_top_proxies/List_of_proteins_for_cox_models_220121.csv")
proteins <- names$Protein %>% as.character()

# Plot each disease by cases and controls vs age 
plot_list <- list()
plot_list2 <- list()

for(i in 1:12){
  disease <- my.list[[i]]
  join <- left_join(disease, data, by = "Sample_Name")
  # remove missing covariates 
  join <- join[!is.na(join$units), ]
  join <- join[!is.na(join$usual), ]
  join <- join[!is.na(join$smokingScore), ]
  join <- join[!is.na(join$simd), ]
  join <- join[!is.na(join$EA), ]
  join <- join[!is.na(join$bmi), ]
  # add status identifier 
  join <- join %>% mutate(Status = case_when(
  Event == "1" ~ "Cases",
  Event == "0" ~ "Controls"))
  # remove NA event rows
  overlap <- which(!join$Event == "NA")
  join <- join[overlap,]
  for (j in proteins){
    id <- j
    dataset <- join
    # get relevant info for plot 
    interest <- dataset[c("Status", id)]
    names(interest)[2] <- "protein"
    # violin plot 
    p <- ggplot(interest, aes(x=Status, y=protein)) +
    geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
    geom_boxplot(width=0.1) + theme_minimal() + ylab(j) + ggtitle(j) + xlab(names[[i]])
    plot_list2[[j]] <- p
  }
  plot_list[[i]] <- plot_list2
}

plot_list[[2]][[1]] # disease then protein


pdf(file = paste0("/Volumes/marioni-lab/Protein_DNAm_Proxies/Plots_for_cox_model_proxy_distributions/All_disease_cases_and_controls_by_proxy_fully_adjusted_model.pdf"))
for (i in 1:12){
  for (j in 1:56){
    print(plot_list[[i]][[j]])
  }
}
dev.off()


### BASIC MODEL 


data <- d1 

my.list = list(AD, bowel, stroke, RA, pain, lung, IHD, diab, depr, IBD, COPD, breast)
names <- c("AD", "bowel", "stroke", "RA", "pain", "lung", "IHD", "diab", "depr", "IBD", "COPD", "breast")

proteins <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Glmnet_predictor_results/Model_results_tables_for_full_training/names_of_56_proteins.csv")
proteins <- proteins$Protein %>% as.character()

# Plot each disease by cases and controls vs age 
plot_list <- list()
plot_list2 <- list()

for(i in 1:12){
  disease <- my.list[[i]]
  join <- left_join(disease, data, by = "Sample_Name")
  # add status identifier 
  join <- join %>% mutate(Status = case_when(
  Event == "1" ~ "Cases",
  Event == "0" ~ "Controls"))
  # remove NA event rows
  overlap <- which(!join$Event == "NA")
  join <- join[overlap,]
  for (j in proteins){
    id <- j
    dataset <- join
    # get relevant info for plot 
    interest <- dataset[c("Status", id)]
    names(interest)[2] <- "protein"
    # violin plot 
    p <- ggplot(interest, aes(x=Status, y=protein)) +
    geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
    geom_boxplot(width=0.1) + theme_minimal() + ylab(j) + ggtitle(j) + xlab(names[[i]])
    plot_list2[[j]] <- p
  }
  plot_list[[i]] <- plot_list2
}

plot_list[[2]][[1]] # disease then protein


pdf(file = paste0("/Volumes/marioni-lab/Protein_DNAm_Proxies/Plots_for_cox_model_proxy_distributions/All_disease_cases_and_controls_by_proxy_basic_model.pdf"))
for (i in 1:12){
  for (j in 1:56){
    print(plot_list[[i]][[j]])
  }
}
dev.off()


### WBC plots - by cases and controls for each disease 


### FULLY ADJUSTED 

data <- d1 

my.list = list(AD, bowel, stroke, RA, pain, lung, IHD, diab, depr, IBD, COPD, breast)
names <- c("AD", "bowel", "stroke", "RA", "pain", "lung", "IHD", "diab", "depr", "IBD", "COPD", "breast")

WBCs <- c("CD4T", "CD8T", "NK", "Bcell", "Gran", "Mono")

# Plot each disease by cases and controls vs age 
plot_list <- list()
plot_list2 <- list()

for(i in 1:12){
  disease <- my.list[[i]]
  join <- left_join(disease, data, by = "Sample_Name")
  # remove missing covariates 
  join <- join[!is.na(join$units), ]
  join <- join[!is.na(join$usual), ]
  join <- join[!is.na(join$smokingScore), ]
  join <- join[!is.na(join$simd), ]
  join <- join[!is.na(join$EA), ]
  join <- join[!is.na(join$bmi), ]
  # add status identifier 
  join <- join %>% mutate(Status = case_when(
  Event == "1" ~ "Cases",
  Event == "0" ~ "Controls"))
  # remove NA event rows
  overlap <- which(!join$Event == "NA")
  join <- join[overlap,]
  for (j in WBCs){
    id <- j
    dataset <- join
    # get relevant info for plot 
    interest <- dataset[c("Status", id)]
    names(interest)[2] <- "WBC"
    # violin plot 
    p <- ggplot(interest, aes(x=Status, y=WBC)) +
    geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
    geom_boxplot(width=0.1) + theme_minimal() + ylab(j) + ggtitle(j) + xlab(names[[i]])
    plot_list2[[j]] <- p
  }
  plot_list[[i]] <- plot_list2
}

plot_list[[2]][[1]] # disease then protein


pdf(file = paste0("/Volumes/marioni-lab/Protein_DNAm_Proxies/Plots_for_cox_model_proxy_distributions/All_disease_cases_and_controls_by_WBC_props_fully_adjusted_model.pdf"))
for (i in 1:12){
  for (j in 1:6){
    print(plot_list[[i]][[j]])
  }
}
dev.off()



### BASIC MODEL

data <- d1 

my.list = list(AD, bowel, stroke, RA, pain, lung, IHD, diab, depr, IBD, COPD, breast)
names <- c("AD", "bowel", "stroke", "RA", "pain", "lung", "IHD", "diab", "depr", "IBD", "COPD", "breast")

WBCs <- c("CD4T", "CD8T", "NK", "Bcell", "Gran", "Mono")

# Plot each disease by cases and controls vs age 
plot_list <- list()
plot_list2 <- list()

for(i in 1:12){
  disease <- my.list[[i]]
  join <- left_join(disease, data, by = "Sample_Name")
  # add status identifier 
  join <- join %>% mutate(Status = case_when(
  Event == "1" ~ "Cases",
  Event == "0" ~ "Controls"))
  # remove NA event rows
  overlap <- which(!join$Event == "NA")
  join <- join[overlap,]
  for (j in WBCs){
    id <- j
    dataset <- join
    # get relevant info for plot 
    interest <- dataset[c("Status", id)]
    names(interest)[2] <- "WBC"
    # violin plot 
    p <- ggplot(interest, aes(x=Status, y=WBC)) +
    geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
    geom_boxplot(width=0.1) + theme_minimal() + ylab(j) + ggtitle(j) + xlab(names[[i]])
    plot_list2[[j]] <- p
  }
  plot_list[[i]] <- plot_list2
}

plot_list[[2]][[1]] # disease then protein


pdf(file = paste0("/Volumes/marioni-lab/Protein_DNAm_Proxies/Plots_for_cox_model_proxy_distributions/All_disease_cases_and_controls_by_WBC_props_basic_model.pdf"))
for (i in 1:12){
  for (j in 1:6){
    print(plot_list[[i]][[j]])
  }
}
dev.off()































