Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# GROUPS AND FOLD ALLOCATIONS ##############################################
############################################################################################
############################################################################################

# Here we assign test/train split and full training individuals into groups based on folds split 
# by plate of DNAm processing batch in the LBC1936 cohort, to account for batch effects 
# This is done for neuro and inflam olink panels 
# Folds files are matched to the order of y variable input to elastic net models 

library(tidyverse)

############################################################################################

### NEUROLOGY FULL TRAIN

############################################################################################

# Load IDs for people that will need to be removed 
ID_remove <- c("LBC360558", "LBC360760", "LBC360536", "LBC361272",
"LBC360213", "LBC360262", "LBC361264", "LBC361030", "LBC360412",
"LBC361076", "LBC360721")

# Read in proteins # 751 people
proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/protein_data/W2_Normalised_Neurology_Proteins.csv")
names(proteins)[1] <- "lbc36no"

# Remove IDs from proteins # 741 people remaining 
overlap <- which(proteins$lbc36no %in% ID_remove)
proteins <- proteins[-overlap,]

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# Filter to wave 2
target <- target %>% filter(WAVE == "2") # 801 people

# Filter to LBC36
target <- target %>% filter(cohort == "LBC36") # 801 people

# Get LBC ID and basenames info 
target <- target[c(2,4,16,18)]
names(target)[2] <- "lbc36no"

# Make sure no factor level issues when joining 
target$lbc36no <- as.character(target$lbc36no)
proteins$lbc36no <- as.character(proteins$lbc36no)

# Work out which IDs are not in the target list for wave 2 
#non_overlap <- which(!proteins$lbc36no %in% target$lbc36no) 

#non_overlap_set <- proteins[non_overlap,]

# Those that dont match up - rob confirmed that these do not have methylation 
#> non_overlap_set$lbc36no
# [1] "LBC360030" "LBC360043" "LBC360087" "LBC360114" "LBC360120" "LBC360128"
# [7] "LBC360161" "LBC360191" "LBC360192" "LBC360194" "LBC360199" "LBC360234"
#[13] "LBC360243" "LBC360275" "LBC360301" "LBC360306" "LBC360311" "LBC360316"
#[19] "LBC360325" "LBC360405" "LBC360420" "LBC360424" "LBC360434" "LBC360439"
#[25] "LBC360465" "LBC360511" "LBC360607" "LBC360659" "LBC360666" "LBC360667"
#[31] "LBC360718" "LBC360754" "LBC360771" "LBC360779" "LBC360790"

# Find which joining IDs in proteins 
overlap <- which(proteins$lbc36no %in% target$lbc36no) 

proteins <- proteins[overlap,] 

# Merge based on only the proteins we have info for
protein <- left_join(proteins, target, by = "lbc36no")

# Reorder for ease 
protein_p <- protein[c(94,95,96,1,2:93)]

# Get a file for ROb for next steps 
# IDS <- protein_p[c(1,4)]

# write to csv for rob 
# write.csv(IDS, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/neuro_analysis/neuro_train_IDS.csv")

# Bin LBC IDs now and keep final protiens 
protein_p <- protein_p[-4]

protein <- protein_p[c(-2,-3)]

protein$Basename <- as.character(protein$Basename)

# Get unique plate names 
plate <- protein_p %>% select(plate) %>% unique()


############################################################################################

### SORT OUT 12 FOLDS AS SEPARATE PLATES FOR CROSS VALIDATION 

### TAKE A LOOK AT NUMBERS PER FOLD FOR ALLOCATIONS

# Dataset that includes the plate info 
plates <- protein_p

# Get just the plate and basenames info 
plates <- plates[c(1,3)]

# Figure out the remaining plates 
remain <- plates %>% select(plate) %>% unique()

### 706 total - to be split into sets 

## Test set = 130
#g2 <- protein_p %>% select(Basename, plate) %>% filter(plate == "6") # 49
#g19 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate22") # 81

## Train set = 576 

# fold 1 (57)
#g5 <- plates %>% select(Basename, plate) %>% filter(plate == "2") #  31
#g6 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate09") # 12
#g7 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate11") #  13
#g9 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate06") # 1

# fold 2 (50)
#g3 <- plates %>% select(Basename, plate) %>% filter(plate == "E141434_LBC-4Meth450_Plate06") # 5
#g4 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate14") # 9
#g8 <- plates %>% select(Basename, plate) %>% filter(plate == "1") # 35
#g30 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate26b") # 1

# fold 3 (59)
#g15 <- plates %>% select(Basename, plate) %>% filter(plate == "9") # 52
#g10 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate08") # 7

# fold 4 (56)
#g1 <- protein_p %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate07") # 12
#g14 <- plates %>% select(Basename, plate) %>% filter(plate == "4") # 44

# fold 5 (58)
#g11 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate17") # 7
#g17 <- plates %>% select(Basename, plate) %>% filter(plate == "8") # 47
#g23 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate13") # 4

# fold 6 (58)
#g16 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate10") # 7
#g20 <- plates %>% select(Basename, plate) %>% filter(plate == "5") # 42
#g18 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate12") # 9

# fold 7 (59)
#g13 <- plates %>% select(Basename, plate) %>% filter(plate == "10") # 43
#g29 <- plates %>% select(Basename, plate) %>% filter(plate == "11") # 12
#g26 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate19") # 3
#g28 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate20") # 1

# fold 8 (68)
#g25 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate21") # 68

# fold 9 (56)
#g27 <- plates %>% select(Basename, plate) %>% filter(plate == "7") # 46
#g21 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate18") # 4
#g24 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate16") # 6

# fold 10 (55)
#g12 <- plates %>% select(Basename, plate) %>% filter(plate == "3") # 44
#g22 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate15") # 11


############################################################################################

### SORT OUT 12 FOLDS

# Dataset that includes the plate info 
plates <- protein_p

# Get just the plate and basenames info 
plates <- plates[c(1,3)]

# Convert to character just in case
plates$plate <- as.character(plates$plate)

# Assign fold ids 
plates <- plates %>% mutate(folds = case_when(
  plate == "2" ~ "1",
  plate == "E11970_Meth450K_Plate09" ~ "1",
  plate == "E11970_Meth450K_Plate11" ~ "1",
  plate == "E11970_Meth450K_Plate06" ~ "1",
  plate == "E141434_LBC-4Meth450_Plate06" ~ "2",
  plate == "E11970_Meth450K_Plate14" ~ "2",
  plate == "1" ~ "2",
  plate == "E11970_Meth450K_Plate26b" ~ "2",
  plate == "9" ~ "3",
  plate == "E11970_Meth450K_Plate08" ~ "3",
  plate == "4" ~ "4",
  plate == "E11970_Meth450K_Plate07" ~ "4",
  plate == "E11970_Meth450K_Plate17" ~ "5",
  plate == "8" ~ "5",
  plate == "E11970_Meth450K_Plate13" ~ "5",
  plate == "E11970_Meth450K_Plate10" ~ "6",
  plate == "5" ~ "6",
  plate == "E11970_Meth450K_Plate12" ~ "6",
  plate == "10" ~ "7",
  plate == "11" ~ "7",
  plate == "E11970_Meth450K_Plate19" ~ "7",
  plate == "E11970_Meth450K_Plate20" ~ "7",
  plate == "E11970_Meth450K_Plate21" ~ "8",
  plate == "7" ~ "9",
  plate == "E11970_Meth450K_Plate18" ~ "9",
  plate == "E11970_Meth450K_Plate16" ~ "9",
  plate == "3" ~ "10",
  plate == "E11970_Meth450K_Plate15" ~ "10",
  plate == "6" ~ "11",
  plate == "E11970_Meth450K_Plate22" ~ "12"))


# Change to character for identical check 
plates$Basename <- as.character(plates$Basename)
protein$Basename <- as.character(protein$Basename)

# Read in y variable 
protein <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train.csv")

# Check that the folds order matches up to the y variable order 
identical(protein$Basename, plates$Basename) # TRUE

# Write out plate folds 
write.csv(plates, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train_folds.csv", row.names = F)

############################################################################################

### NEURO TEST TRAIN SPLIT 

############################################################################################

# Load IDs for people that will need to be removed 
ID_remove <- c("LBC360558", "LBC360760", "LBC360536", "LBC361272",
"LBC360213", "LBC360262", "LBC361264", "LBC361030", "LBC360412",
"LBC361076", "LBC360721")

# Read in proteins # 751 people
proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/protein_data/W2_Normalised_Neurology_Proteins.csv")
names(proteins)[1] <- "lbc36no"

# Remove IDs from proteins # 741 people remaining 
overlap <- which(proteins$lbc36no %in% ID_remove)
proteins <- proteins[-overlap,]

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# Filter to wave 1 
target <- target %>% filter(WAVE == "2") # 801 people

# Filter to LBC36
target <- target %>% filter(cohort == "LBC36") # 801 people

# Get LBC ID and basenames info 
target <- target[c(2,4,16,18)]
names(target)[2] <- "lbc36no"

# Make sure no factor level issues when joining 
target$lbc36no <- as.character(target$lbc36no)
proteins$lbc36no <- as.character(proteins$lbc36no)

# Work out which IDs are not in the target list for wave 2 
#non_overlap <- which(!proteins$lbc36no %in% target$lbc36no) 

#non_overlap_set <- proteins[non_overlap,]

# Those that dont match up - rob confirmed that these do not have methylation 
#> non_overlap_set$lbc36no
# [1] "LBC360030" "LBC360043" "LBC360087" "LBC360114" "LBC360120" "LBC360128"
# [7] "LBC360161" "LBC360191" "LBC360192" "LBC360194" "LBC360199" "LBC360234"
#[13] "LBC360243" "LBC360275" "LBC360301" "LBC360306" "LBC360311" "LBC360316"
#[19] "LBC360325" "LBC360405" "LBC360420" "LBC360424" "LBC360434" "LBC360439"
#[25] "LBC360465" "LBC360511" "LBC360607" "LBC360659" "LBC360666" "LBC360667"
#[31] "LBC360718" "LBC360754" "LBC360771" "LBC360779" "LBC360790"

# Find which joining IDs in proteins 
overlap <- which(proteins$lbc36no %in% target$lbc36no) 

proteins <- proteins[overlap,] 

# Merge based on only the proteins we have info for
protein <- left_join(proteins, target, by = "lbc36no")

# Reorder for ease 
protein_p <- protein[c(94,95,96,1,2:93)]

# Bin LBC IDs now and keep final protiens 
protein_p <- protein_p[-4]

protein <- protein_p[c(-2,-3)]

protein$Basename <- as.character(protein$Basename)

# Get unique plate names 
plate <- protein_p %>% select(plate) %>% unique()



############################################################################################

### SORT OUT 12 FOLDS AS SEPARATE PLATES FOR CROSS VALIDATION 

############################################################################################

### TAKE A LOOK AT NUMBERS PER FOLD FOR ALLOCATIONS

# Dataset that includes the plate info 
plates <- protein_p

# Get just the plate and basenames info 
plates <- plates[c(1,3)]

# Figure out the remaining plates 
remain <- plates %>% select(plate) %>% unique()


#                           plate
#1        E11970_Meth450K_Plate07
#2                              6
#7   E141434_LBC-4Meth450_Plate06
#10       E11970_Meth450K_Plate14
#11                             2
#12       E11970_Meth450K_Plate09
#13       E11970_Meth450K_Plate11
#14                             1
#17       E11970_Meth450K_Plate06
#23       E11970_Meth450K_Plate08
#28       E11970_Meth450K_Plate17
#29                             3
#40                            10
#41                             4
#50                             9
#52       E11970_Meth450K_Plate10
#76                             8
#79       E11970_Meth450K_Plate12
#89       E11970_Meth450K_Plate22
#95                             5
#96       E11970_Meth450K_Plate18
#99       E11970_Meth450K_Plate15
#103      E11970_Meth450K_Plate13
#139      E11970_Meth450K_Plate16
#150      E11970_Meth450K_Plate21
#162      E11970_Meth450K_Plate19
#185                            7
#206      E11970_Meth450K_Plate20
#420                           11
#474     E11970_Meth450K_Plate26b


### 706 total - to be split into sets 

## Test set = 130
g2 <- protein_p %>% select(Basename, plate) %>% filter(plate == "6") # 49
g19 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate22") # 81

## Train set = 576 

# fold 1 (57)
g5 <- plates %>% select(Basename, plate) %>% filter(plate == "2") #  31
g6 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate09") # 12
g7 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate11") #  13
g9 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate06") # 1

# fold 2 (50)
g3 <- plates %>% select(Basename, plate) %>% filter(plate == "E141434_LBC-4Meth450_Plate06") # 5
g4 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate14") # 9
g8 <- plates %>% select(Basename, plate) %>% filter(plate == "1") # 35
g30 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate26b") # 1

# fold 3 (59)
g15 <- plates %>% select(Basename, plate) %>% filter(plate == "9") # 52
g10 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate08") # 7

# fold 4 (56)
g1 <- protein_p %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate07") # 12
g14 <- plates %>% select(Basename, plate) %>% filter(plate == "4") # 44

# fold 5 (58)
g11 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate17") # 7
g17 <- plates %>% select(Basename, plate) %>% filter(plate == "8") # 47
g23 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate13") # 4

# fold 6 (58)
g16 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate10") # 7
g20 <- plates %>% select(Basename, plate) %>% filter(plate == "5") # 42
g18 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate12") # 9

# fold 7 (59)
g13 <- plates %>% select(Basename, plate) %>% filter(plate == "10") # 43
g29 <- plates %>% select(Basename, plate) %>% filter(plate == "11") # 12
g26 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate19") # 3
g28 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate20") # 1

# fold 8 (68)
g25 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate21") # 68

# fold 9 (56)
g27 <- plates %>% select(Basename, plate) %>% filter(plate == "7") # 46
g21 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate18") # 4
g24 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate16") # 6

# fold 10 (55)
g12 <- plates %>% select(Basename, plate) %>% filter(plate == "3") # 44
g22 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate15") # 11


############################################################################################

### SORT OUT TEST SET

############################################################################################

# Test set 
test_group_1 <- protein_p %>% select(Basename, plate) %>% filter(plate == "6") # 49
test_group_2 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate22") # 81
join_group <- rbind(test_group_1, test_group_2) # 130
list1 <- join_group$Basename %>% as.data.frame() # 130
names(list1)[1] <- "group"

############################################################################################

### SORT OUT FOLDS FOR TRAINING SET

############################################################################################

# Dataset that includes the plate info 
plates <- protein_p

# Remove the test set from the data so we just have train data 
overlap <- which(plates$Basename %in% list1$group)
plates <- plates[-overlap,]

# Get just the plate and basenames info 
plates <- plates[c(1,3)]

# Convert to character just in case
plates$plate <- as.character(plates$plate)

# Assign fold ids 
plates <- plates %>% mutate(folds = case_when(
  plate == "2" ~ "1",
  plate == "E11970_Meth450K_Plate09" ~ "1",
  plate == "E11970_Meth450K_Plate11" ~ "1",
  plate == "E11970_Meth450K_Plate06" ~ "1",
  plate == "E141434_LBC-4Meth450_Plate06" ~ "2",
  plate == "E11970_Meth450K_Plate14" ~ "2",
  plate == "1" ~ "2",
  plate == "E11970_Meth450K_Plate26b" ~ "2",
  plate == "9" ~ "3",
  plate == "E11970_Meth450K_Plate08" ~ "3",
  plate == "4" ~ "4",
  plate == "E11970_Meth450K_Plate07" ~ "4",
  plate == "E11970_Meth450K_Plate17" ~ "5",
  plate == "8" ~ "5",
  plate == "E11970_Meth450K_Plate13" ~ "5",
  plate == "E11970_Meth450K_Plate10" ~ "6",
  plate == "5" ~ "6",
  plate == "E11970_Meth450K_Plate12" ~ "6",
  plate == "10" ~ "7",
  plate == "11" ~ "7",
  plate == "E11970_Meth450K_Plate19" ~ "7",
  plate == "E11970_Meth450K_Plate20" ~ "7",
  plate == "E11970_Meth450K_Plate21" ~ "8",
  plate == "7" ~ "9",
  plate == "E11970_Meth450K_Plate18" ~ "9",
  plate == "E11970_Meth450K_Plate16" ~ "9",
  plate == "3" ~ "10",
  plate == "E11970_Meth450K_Plate15" ~ "10"))


# Change to character for identical check 
plates$Basename <- as.character(plates$Basename)

# Read in y variable 
y <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_train_only.csv")

# Check that the folds order matches up to the y variable order 
identical(y$Basename, plates$Basename) # TRUE

# Write out plates 
write.csv(plates, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_train_only_folds.csv", row.names = F)


############################################################################################

### INFLAM FULL TRAIN

############################################################################################


# Load IDs for people that will need to be removed 
ID_remove <- c("LBC360558", "LBC360760", "LBC360536", "LBC361272",
"LBC360213", "LBC360262", "LBC361264", "LBC361030", "LBC360412",
"LBC361076", "LBC360721")

# Read in 70 proteins 
proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/protein_data/Normalised_Inflammatory_Proteins_W1.csv")

# Remove plate ID (not needed as feature)
proteins <- proteins[-72]

# Remove IDs from proteins 
overlap <- which(proteins$lbc36no %in% ID_remove)
proteins <- proteins[-overlap,]

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# Filter to wave 1 
target <- target %>% filter(WAVE == "1") # 1342

# Filter to LBC36
target <- target %>% filter(cohort == "LBC36") # 906

# Get LBC ID and basenames info 
target <- target[c(2,4,16,18)]
names(target)[2] <- "lbc36no"

# Make sure no factor level issues when joining 
target$lbc36no <- as.character(target$lbc36no)
proteins$lbc36no <- as.character(proteins$lbc36no)

# Find which joining IDs in proteins 
overlap <- which(proteins$lbc36no %in% target$lbc36no) 

# Subset target file by these names
proteins <- proteins[overlap,] 

# Merge based on only the proteins we have info for
protein <- left_join(proteins, target, by = "lbc36no")

# Reorder for ease 
protein_p <- protein[c(72,73,74,1,2:71)]

# Bin LBC IDs now and keep final protiens 
protein_p <- protein_p[-4]

protein <- protein_p[c(-2,-3)]

protein$Basename <- as.character(protein$Basename)

# Get unique plate names 
plate <- protein_p %>% select(plate) %>% unique()


############################################################################################

### SORT OUT 12 FOLDS AS SEPARATE PLATES FOR CROSS VALIDATION 

############################################################################################

### TAKE A LOOK AT NUMBERS PER FOLD FOR ALLOCATIONS

# Dataset that includes the plate info 
plates <- protein_p

# Get just the plate and basenames info 
plates <- plates[c(1,3)]

# Figure out the remaining plates 
remain <- plates %>% select(plate) %>% unique()

# These were the test/train previously but can now be added as extra 2 folds 
test_group_1 <- protein_p %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate17") # 72 - now 71
test_group_2 <- protein_p %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate16") # 80 - now 79 

# These I can make 9 folds with 
#p10_10 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate10") # 75
#p11_11 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate11") # 71 - now 69
#p13 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate13") # 63 - now 62 
#p12 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate12") # 73 - now 73
#p14 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate14") # 71 - now 69 
#p15 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate15") # 75 - now 75
#p18 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate18") # 86 - now 85
#p19 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate19") # 71 - now 70
#p20 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate20") # 66 - now 66

# These need to be merged to an extra fold (adds up to 83) 
#p8 <- plates %>% select(Basename, plate) %>% filter(plate == "8") # 5
#p10 <- plates %>% select(Basename, plate) %>% filter(plate == "10") # 4
#p11 <- plates %>% select(Basename, plate) %>% filter(plate == "11") # 3
#p6 <- plates %>% select(Basename, plate) %>% filter(plate == "6") # 4
#p2 <- plates %>% select(Basename, plate) %>% filter(plate == "2") # 16
#p4 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate04") # 8
#p3 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate03") # 8 - now 7
#p06 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate06") # 15
#p21 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate21") # 11 - now 10
#p05 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate05") # 9

# Convert to character just in case
plates$plate <- as.character(plates$plate)

# Assign fold ids 
plates <- plates %>% mutate(folds = case_when(
  plate == "E11970_Meth450K_Plate10" ~ "1",
  plate == "E11970_Meth450K_Plate11" ~ "2",
  plate == "E11970_Meth450K_Plate13" ~ "3",
  plate == "E11970_Meth450K_Plate12" ~ "4",
  plate == "E11970_Meth450K_Plate14" ~ "5",
  plate == "E11970_Meth450K_Plate15" ~ "6",
  plate == "E11970_Meth450K_Plate18" ~ "7",
  plate == "E11970_Meth450K_Plate19" ~ "8",
  plate == "E11970_Meth450K_Plate20" ~ "10",
  plate == "E11970_Meth450K_Plate04" ~ "9",
  plate == "E11970_Meth450K_Plate03" ~ "9",
  plate == "E11970_Meth450K_Plate06" ~ "9",
  plate == "E11970_Meth450K_Plate21" ~ "9",
  plate == "E11970_Meth450K_Plate05" ~ "9",
  plate == "8" ~ "9",
  plate == "10" ~ "9",
  plate == "11" ~ "9",
  plate == "6" ~ "9",
  plate == "2" ~ "9",
  plate == "E11970_Meth450K_Plate17" ~ "11",
  plate == "E11970_Meth450K_Plate16" ~ "12"))

# Change to character for identical check 
plates$Basename <- as.character(plates$Basename)

# Check that the folds order matches up to the y variable order 
identical(protein$Basename, plates$Basename) # TRUE

# Restrict so you have the folds for the model provided as an integer string of values to input into function
folds <- plates$folds %>% as.data.frame()
names(folds)[1] <- "F"
folds$F <- as.integer(folds$F)
test <- folds$F # 875 long now


# Change to character for identical check 
plates$Basename <- as.character(plates$Basename)

# Read in y variable 
y <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train.csv")

# Check that the folds order matches up to the y variable order 
identical(y$Basename, plates$Basename) # TRUE

# Write out plates 
write.csv(plates, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train_folds.csv", row.names = F)



############################################################################################

### INFLAM TEST TRAIN SPLIT 

############################################################################################

# Load IDs for people that will need to be removed 
ID_remove <- c("LBC360558", "LBC360760", "LBC360536", "LBC361272",
"LBC360213", "LBC360262", "LBC361264", "LBC361030", "LBC360412",
"LBC361076", "LBC360721")

# Read in 70 proteins 
proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/protein_data/Normalised_Inflammatory_Proteins_W1.csv")

# Remove plate ID (not needed as feature)
proteins <- proteins[-72]

# Remove IDs from proteins 
overlap <- which(proteins$lbc36no %in% ID_remove)
proteins <- proteins[-overlap,]

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# Filter to wave 1 
target <- target %>% filter(WAVE == "1") # 1342

# Filter to LBC36
target <- target %>% filter(cohort == "LBC36") # 906

# Get LBC ID and basenames info 
target <- target[c(2,4,16,18)]
names(target)[2] <- "lbc36no"

# Make sure no factor level issues when joining 
target$lbc36no <- as.character(target$lbc36no)
proteins$lbc36no <- as.character(proteins$lbc36no)

# Find which joining IDs in proteins 
overlap <- which(proteins$lbc36no %in% target$lbc36no) 

# Subset target file by these names
proteins <- proteins[overlap,] 

# Merge based on only the proteins we have info for
protein <- left_join(proteins, target, by = "lbc36no")

# Reorder for ease 
protein_p <- protein[c(72,73,74,1,2:71)]


# Bin LBC IDs now and keep final protiens 
protein_p <- protein_p[-4]

protein <- protein_p[c(-2,-3)]

protein$Basename <- as.character(protein$Basename)

# Get unique plate names 
plate <- protein_p %>% select(plate) %>% unique()





############################################################################################

### SORT OUT TEST SET OF 152 PROTEINS FROM 2 PLATES 

############################################################################################

# Okay so 71 and 79 = 150 which i can use as first test set for now (plates 17 and 16)
test_group_1 <- protein_p %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate17") # 71
test_group_2 <- protein_p %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate16") # 79
join_group <- rbind(test_group_1, test_group_2) # 152 
list1 <- join_group$Basename %>% as.data.frame() # 152
names(list1)[1] <- "group"


# Now I need to subset my ytrain and ytest based on these 2 groups
y1 <- which(protein$Basename %in% list1$group)
y1 <- protein[y1,] # 150
y1$Basename <- as.character(y1$Basename)
y2 <- which(!protein$Basename %in% list1$group)
y2 <- protein[y2,] # 725
y2$Basename <- as.character(y2$Basename)



############################################################################################

### SORT OUT TRAINING SET (734) GROUPED BY PLATES TO AVOID LEAKAGE 

############################################################################################

### TAKE A LOOK AT REMAINING FOLDS FOR TRAINING SET 

# Dataset that includes the plate info 
plates <- protein_p

# Remove the test set from the data so we just have train data 
overlap <- which(plates$Basename %in% list1$group)
plates <- plates[-overlap,]

# Get just the plate and basenames info 
plates <- plates[c(1,3)]

# Figure out the remaining plates 
remain <- plates %>% select(plate) %>% unique()

# These I can make 9 folds with 
p10_10 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate10") # 75
p11_11 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate11") # 71
p13 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate13") # 63
p12 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate12") # 73
p14 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate14") # 71
p15 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate15") # 75
p18 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate18") # 86
p19 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate19") # 71
p20 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate20") # 66

# These need to be merged to an extra fold (adds up to 83) 
p8 <- plates %>% select(Basename, plate) %>% filter(plate == "8") # 5
p10 <- plates %>% select(Basename, plate) %>% filter(plate == "10") # 4
p11 <- plates %>% select(Basename, plate) %>% filter(plate == "11") # 3
p6 <- plates %>% select(Basename, plate) %>% filter(plate == "6") # 4
p2 <- plates %>% select(Basename, plate) %>% filter(plate == "2") # 16
p4 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate04") # 8
p3 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate03") # 8
p06 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate06") # 15
p21 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate21") # 11
p05 <- plates %>% select(Basename, plate) %>% filter(plate == "E11970_Meth450K_Plate05") # 9

# Convert to character just in case
plates$plate <- as.character(plates$plate)

# Assign fold ids 
plates <- plates %>% mutate(folds = case_when(
  plate == "E11970_Meth450K_Plate10" ~ "1",
  plate == "E11970_Meth450K_Plate11" ~ "2",
  plate == "E11970_Meth450K_Plate13" ~ "3",
  plate == "E11970_Meth450K_Plate12" ~ "4",
  plate == "E11970_Meth450K_Plate14" ~ "5",
  plate == "E11970_Meth450K_Plate15" ~ "6",
  plate == "E11970_Meth450K_Plate18" ~ "7",
  plate == "E11970_Meth450K_Plate19" ~ "8",
  plate == "E11970_Meth450K_Plate20" ~ "10",
  plate == "E11970_Meth450K_Plate04" ~ "9",
  plate == "E11970_Meth450K_Plate03" ~ "9",
  plate == "E11970_Meth450K_Plate06" ~ "9",
  plate == "E11970_Meth450K_Plate21" ~ "9",
  plate == "E11970_Meth450K_Plate05" ~ "9",
  plate == "8" ~ "9",
  plate == "10" ~ "9",
  plate == "11" ~ "9",
  plate == "6" ~ "9",
  plate == "2" ~ "9"))

# Check that the folds order matches up to the ytrain variable order 
plates$Basename <- as.character(plates$Basename)
identical(y2$Basename, plates$Basename) # TRUE

# Restrict so you have the folds for the model provided as an integer string of values to input into function
folds <- plates$folds %>% as.data.frame()
names(folds)[1] <- "F"
folds$F <- as.integer(folds$F)
test <- folds$F


# Change to character for identical check 
plates$Basename <- as.character(plates$Basename)

# Read in y variable 
y <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_train_only.csv")

# Check that the folds order matches up to the y variable order 
identical(y$Basename, plates$Basename) # TRUE

# Write out plates 
write.csv(plates, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_train_only_folds.csv", row.names = F)



############################################################################################
