Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.


##################################################################################################

### MAKE THE FORREST PLOT FOR NON DIABETES / COPD 

##################################################################################################

# setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

library(tidyverse)


x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")

names1 <- unique(x$Name)

# create a key that can be used for the plot naming - removing nTRK3 and CLEC11A repeats (i.e. 78 unique predictor names)
y <- x[c(1,4,25,26)] %>% as.data.frame()

# save this out to send to rob for the suppl figures so we are consistent 
# write.csv(y, "Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv", row.names = F)


##################################################################

# Plot the forest for each of the disease outcomes and then join using patchwork 

library(tidyverse)
library(ggplot2)
library(readxl)
library(patchwork)

setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

library(tidyverse)

# load in the datset withthe HR info for FA results 
x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")
x <- as.data.frame(x)

# Join in the correct naming for the plots as generated above 
naming <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
names(naming)[4] <- "Plot_names"
naming <- naming[c(1,4)]
naming$Predictor <- as.character(naming$Predictor)

# check order of predictor is the same 
identical(naming$Predictor, x$Predictor) # TRUE 

# Join in the plot naming 
naming <- naming[2]
data <- cbind(x, naming)

### Pain
x <- data
x <- x[which(x$Outcome %in% c("Pain")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.title.x = element_text(size = 27),
  axis.text.x = element_text(size = 27),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

pain <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme


pdf("pdf_pain.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()


### RA
x <- data
x <- x[which(x$Outcome %in% c("RA")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

RA <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme


pdf("pdf_RA.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()


### Dep
x <- data
x <- x[which(x$Outcome %in% c("Depression")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.title.x = element_text(size = 27),
  axis.text.x = element_text(size = 27),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

Dep <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme


pdf("pdf_Dep.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()


### Bowel
x <- data
x <- x[which(x$Outcome %in% c("Bowel.Cancer")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

x$Outcome2 <- "Bowel Cancer"

My_Theme = theme(
  axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")


pdf("pdf_Bowel.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome2, scales = "free_y") + My_Theme
  dev.off()

### AD
x <- data
x <- x[which(x$Outcome %in% c("Alzheimer's Disease")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

x$Outcome2 <- "Alzheimer's Dementia"

My_Theme = theme(
  axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

AD <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome2, scales = "free_y") + My_Theme

pdf("pdf_AD.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome2, scales = "free_y") + My_Theme
  dev.off()


### Lung
x <- data
x <- x[which(x$Outcome %in% c("Lung.Cancer")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

x$Outcome2 <- "Lung Cancer"

My_Theme = theme(
  axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

Lung <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome2, scales = "free_y") + My_Theme


pdf("pdf_Lung.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome2, scales = "free_y") + My_Theme
  dev.off()



### IHD
x <- data
x <- x[which(x$Outcome %in% c("IHD")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

IHD <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

pdf("pdf_IHD.pdf", width = 7, height = 7)
 ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()



### IBD
x <- data
x <- x[which(x$Outcome %in% c("IBD")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

IBD <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

pdf("pdf_IBD.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()


### Stroke
x <- data
x <- x[which(x$Outcome %in% c("Stroke")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.title.x = element_text(size = 27),
  axis.text.x = element_text(size = 27),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

Stroke <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

pdf("pdf_Stroke.pdf", width = 7, height = 7)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("")+ xlab ("")+ ylim(0.3, 2.7) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()


### COPD
x <- data
x <- x[which(x$Outcome %in% c("COPD")),]

# Assign x and y variables for forest plot 
x$Outcome2 <- x$Outcome
# levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  axis.title.x = element_text(size = 27),
  axis.text.x = element_text(size = 27),
  axis.text.y = element_text(size = 27),
  axis.title.y = element_text(size = 32),
  strip.text = element_text(size = 30, face = "bold"),
  legend.text=element_text(size=27),
  legend.title=element_text(size=32, face = "bold"), legend.position = "none")

COPD <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+ ylim(0.3, 2.7) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

pdf("pdf_COPD.pdf", width = 10, height = 21)
ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 4.5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+ ylim(0.3, 2.7) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()


# Will have to piece together in inkscape - the patchwork isnt collatign them as per the original!

# # Now try to patchwork them together 

# plot1 <- (AD + Bowel + RA) / (IBD + IHD + Lung) / (pain + Dep + Stroke)

# ### COMBINE PLOTS 1 and 2 
# library(patchwork)
# pdf("testplot3.pdf", width = 35, height = 25)
# plot1 + COPD + 
#   plot_layout(widths = c(2, 1))
# dev.off()




# # Now plot the forest plot for the associations with diabetes removed 

# setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

# library(tidyverse)

# # load in the datset withthe HR info for FA results 
# x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")
# x <- as.data.frame(x)

# # Join in the correct naming for the plots as generated above 
# naming <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
# names(naming)[4] <- "Plot_names"
# naming <- naming[c(1,4)]
# naming$Predictor <- as.character(naming$Predictor)

# # check order of predictor is the same 
# identical(naming$Predictor, x$Predictor) # TRUE 

# # Join in the plot naming 
# naming <- naming[2]
# x <- cbind(x, naming)

# # Remove COPD and diabetes 
# x <- x[-which(x$Outcome %in% c("COPD", "Diabetes")),]


# # Change disease names so they are on 2 lines 
# data <- x

# data <- data %>%
# mutate(Outcome = str_replace(Outcome, "Breast.Cancer", "Breast Cancer"))
# data <- data %>%
# mutate(Outcome = str_replace(Outcome, "Lung.Cancer", "Lung Cancer"))
# data <- data %>%
# mutate(Outcome = str_replace(Outcome, "Bowel.Cancer", "Bowel Cancer"))
# data <- data %>%
# mutate(Outcome = str_replace(Outcome, "Alzheimer's Disease", "Alzheimer's Dementia"))
# data <- data %>%
# mutate(Outcome = str_replace(Outcome, "RA", "Rheumatoid Arthritis"))
# data <- data %>%
# mutate(Outcome = str_replace(Outcome, "IHD", "IHD"))

# x <- data 

# My_Theme = theme(
#   axis.title.x = element_text(size = 27),
#   axis.text.x = element_text(size = 27),
#   axis.text.y = element_text(size = 27),
#   axis.title.y = element_text(size = 32),
#   strip.text = element_text(size = 30, face = "bold"),
#   legend.text=element_text(size=27),
#   legend.title=element_text(size=32, face = "bold"), legend.position = "none")

# # Assign x and y variables for forest plot 
# x$Outcome2 <- x$Outcome
# # levels(x$Outcome2) = c(" ", "  ", "    ")
# x$TraitVar <- paste0(x$Plot_names)
# x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))


# ### SAVE THIS PLOT AS A PLOT1 VARIABLE FOR JOINING TO MAIN COPD 

# # plot1 <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
# #   geom_point(position=position_dodge(width=0.5), shape = 20, size = 7)+
# #   geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
# #                 position = position_dodge(0.5), width = 0.05,
# #                 colour = "black")+
# #   ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
# #   geom_hline(yintercept = 1, linetype = "dotted")+
# #   theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "right",
# #         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
# #   coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

# # test without the theme 
#  ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
#   geom_point(position=position_dodge(width=0.5), shape = 20, size = 7)+
#   geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
#                 position = position_dodge(0.5), width = 0.05,
#                 colour = "black")+
#   ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "right",
#         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
#   coord_flip() + facet_wrap(~Outcome, scales = "free_y")


#  ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
#   geom_point(position=position_dodge(width=0.5), shape = 20, size = 7)+
#   geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
#                 position = position_dodge(0.5), width = 0.05,
#                 colour = "black")+
#   ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "right",
#         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
#   coord_flip() + facet_wrap(~Outcome, scales = "free_y")


# ##################################################################################################################

# ### NOW DO COPD PLOT 
# library(tidyverse)
# library(ggplot2)
# library(readxl)
# library(patchwork)


# setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

# library(tidyverse)

# # load in the datset withthe HR info for FA results 
# x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")
# x <- as.data.frame(x)

# # Join in the correct naming for the plots as generated above 
# naming <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
# names(naming)[4] <- "Plot_names"
# naming <- naming[c(1,4)]
# naming$Predictor <- as.character(naming$Predictor)

# # check order of predictor is the same 
# identical(naming$Predictor, x$Predictor) # TRUE 

# # Join in the plot naming 
# naming <- naming[2]
# x <- cbind(x, naming)

# # Remove COPD and diabetes 
# x <- x[which(x$Outcome %in% c("COPD")),]


# # Assign x and y variables for forest plot 
# x$Outcome2 <- x$Outcome
# # levels(x$Outcome2) = c(" ", "  ", "    ")
# x$TraitVar <- paste0(x$Plot_names)
# x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))


# # order by hazard ratio 
# # x<- x[rev(order(x$Hazard.Ratio.x)),]

# # a <- ifelse(x$Panel == "SomaScan", "red", "blue")

# My_Theme = theme(
#   axis.title.x = element_text(size = 27),
#   axis.text.x = element_text(size = 27),
#   axis.text.y = element_text(size = 27),
#   axis.title.y = element_text(size = 32),
#   strip.text = element_text(size = 30, face = "bold"),
#   legend.text=element_text(size=27),
#   legend.title=element_text(size=32, face = "bold"), legend.position = "none")


# plot2 <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
#   geom_point(size = 4.5)+
#   geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
#                 position = position_dodge(0.5), width = 0.05,
#                 colour = "black")+
#   ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
#         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
#   coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

# ### COMBINE PLOTS 1 and 2 
# library(patchwork)
# pdf("testplot.pdf", width = 35, height = 25)
# plot1 + plot2 + 
#   plot_layout(widths = c(2, 1))
# dev.off()

