####################################################################################

# Covid EpiScores analyses 

####################################################################################

library(tidyverse)
library(readxl)
library(lme4)

# Load covid data from archie
t <- read_excel("C:/Users/s1888864/Desktop/Offline_train/2021-09-03 C19 cases Feb.xlsx")
t <- as.data.frame(t)

# Take a look at variables available in this file (including those without DNAm data)
length(unique(t$id)) # 1713 unique individuals 
table(t$covid) # 554 with diagnosis of covid 
table(t$smr) # 31 in hospital 
table(t$icu) # 6 in ICU 

e <- t %>% filter(t$covid == "1") 
names(e)[1] <- "Sample_Name"

# Load the episcores d1 file that is rank transformed 
d1 <- read.csv("C:/Users/s1888864/Desktop/Offline_train/d1_080321.csv", check.names = F)

# Join episcores file into those with covid 
d2 <- left_join(e, d1, by = "Sample_Name")
dim(d2)

# How many individuals have episcore data 
table(is.na(d2$SMPD1))
# 269 people have covid and episcores data 
table(d2$smr)
# 29 cases


library(glm2)

d2$smr <- as.factor(d2$smr)
results <- data.frame(episcore = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers <- colnames(d2)[117:226]
outcome <- "hospitalised covid"

for (i in 1:length(markers)){
     name <- as.character(markers[i])
     m <- glm(smr ~ scale(d2[,name]) + scale(Age) + factor(Female), data = d2, 
     family = binomial)

     Beta <- coef(summary(m))[2,1]
     SE <- coef(summary(m))[2,2]
     p <- coef(summary(m))[2,4]
     n <- nobs(m)

     results[i,1] <- name
     results[i,2] <- outcome
     results[i,3] <- n
     results[i,4] <- Beta
     results[i,5] <- SE
     results[i,6] <- p

     print(name)
     print(outcome)
}


# Order results by P 
result <- results[order(results$p),]

top <- result[which(result$p < 0.05),] # 6 associations - which dont have the trigger for singular fits 



# Add annotation context 
anno <- read_excel("C:/Users/s1888864/Desktop/Offline_train/Annotations_for_reference.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,4,18,13)]
names(anno)[1] <- "episcore"

result <- left_join(result, anno, by = "episcore")


# Write off results file 
write.csv(result, "C:/Users/s1888864/Desktop/Offline_train/result_glm_covid_only_060921.csv", row.names = F)


## Check family associations present in the dataset 

# Read in the prepped file to cluster 
ped <- read.csv("C:/Users/s1888864/Desktop/Offline_train/pedigree_formatted.csv")

# Join pedigree info to the main dataset as per example above 
names(ped)[2] <- "Sample_Name"
d3 <- left_join(d2, ped, by = "Sample_Name")

d4 <- d3[which(d3$SMPD1 != "NA"),]

length(unique(d4$famid))

table(d4$famid)

t <- d4[which(d4$famid == "1885"),]

## Join both results files for long covid and hospitalisations 

library(tidyverse)
library(readxl)

hosp <- read.csv("Y:/Danni/COVID/result_glm_covid_only_060921.csv")

long <- read.csv("Y:/Danni/COVID/result_glm_long_covid_150921.csv")

# match the order of episcores in the long covid file to the hosp file 

long <- long[match(hosp$episcore, long$episcore),]

# remove uniprot annotations 

hosp <- hosp[c(1:6)]
long <- long[c(1:6)]

# join results together 

res <- cbind(hosp, long)
write.csv(res, "Y:/Danni/COVID/result_suppl_table_covid_joint_131021.csv", row.names = F)


