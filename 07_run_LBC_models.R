Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# Elastic net models for LBC ###############################################
############################################################################################
############################################################################################

# Here, we train elastic net models for Olink LBC1936 protein levels in LBC1936 
# We conduct train-only runs on the full LBC1936 sets for inflam or neuro panels 
# We conduct test/train splits for a holdout assessment of performance

### INFLAM TRAIN FULL

# Required libraries
library(glmnet)
library(tidyverse)
library(foreign)

# Inputs
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train.rds")
ytrain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_full_train.csv")
folds <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train_folds.csv")
output_location <- "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/inflam_full_train/"
cv <- 12

# Sanity check 
identical(folds$Basename, ytrain$Basename)
identical(rownames(xtrain), ytrain$Basename)

# Load model
analysis <- function(xtrain, ytrain, cv, folds) {
 
  set.seed(1783) # Set seed
  seed <- 1783

  Output <- list() # Create storage for predictor weights results
  Output2 <- list() # Create storage for the features table summary 

  for (i in cols) { # Loop through the proteins by columns index which can be set below, before running the analysis function
    tryCatch({
      message("1. Assigning the protein variable.")
      x <- xtrain
      q <- ytrain[1] # Get just basenames for people in the y variable 
      p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
      name_p <- colnames(p) # Get the name of the protein for this iteration
      y <- cbind(q,p) # Bind Basename and protein data together into one set 
      index <- identical(y$Basename, rownames(x)) # check this has been done
      names(y)[2] <- "protein" # Assign a generic name to the protein variable
      y <- as.numeric(y$protein) # Create a numeric list for this variable to feed in as y to the model

      message(paste0("2. Computing predictor weights for ", name_p)) 
      fold <- folds$folds
      lasso.cv <- cv.glmnet(x, y, family="gaussian", alpha = 0.5, foldid = fold, nfolds=cv) # cross validation to get best lambda
      fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, standardize = F, lambda = lasso.cv$lambda.min) # model fit 
      coefs <- coef(fit) # Extract coeficients 
      coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
      coefs$Predictor <- name_p # Assign protein identifier
      names(coefs)[1] <- "Coefficient" # Tidy naming 
      coefs$CpG_Site <- rownames(coefs) # Create cpg column
      coefs <- coefs[c(3,1,2)] # order 
      coefs3 <- coefs[-1,] # Remove intercept (if there was only intercept, this will now have 0 rows)
      write.csv(coefs3, file = paste0(output_location, "01_predictor_weights/predictor_weights_for_protein_", name_p, ".csv"), row.names = F)

      features <- nrow(coefs3) # Count number of features contributing for the protein - i.e. if its > 0 it has features 
      Statement <- matrix(nrow = (length(cols)), ncol = 2) %>% as.data.frame() # Create a vessel to store the results from the 

      if (features > 0) {
        rownames(coefs3) <- NULL # remove rows as they arent needed 
        results <- coefs3 # save weights as outputs list to bind together for all proteins
        message(paste0("3. Your weights for protein ", name_p, " have been calculated."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
        
      } else {
        results <- NULL
        message(paste0("3. Your model for protein ", name_p, " did not generate features."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
      }
      Output[[i]] <- results
      Output2[[i]] <- Statement
    }, error = function(e) cat("skipped"))
  }

  weight_total <- do.call(rbind, Output) # Bind predictor weights together
  write.csv(weight_total, paste0(output_location, "weights_output_210121.csv"), row.names = F)

  statement_total <- do.call(rbind, Output2) # Bind summary table together 
  statement_total <- na.omit(statement_total)
  names(statement_total) <- c("Uniprot_ID", "Features")
  write.csv(statement_total, paste0(output_location, "statement_output_210121.csv"), row.names = F)
  
  return(statement_total)
  message("Model running complete! The outputs will be in your specified location.")
}

message("6. Run the model.")

cols <- c(2:71) # Define cols to loop through here for the protein variables 
cv <- 12 # tell it which cv to run for 

start_time <- Sys.time()
results <- analysis(xtrain = xtrain, ytrain = ytrain, cv = cv, folds = folds)
end_time <- Sys.time()
time_difference <- end_time - start_time 


############################################################################################

### NEURO TRAIN FULL

# Required libraries
library(glmnet)
library(tidyverse)
library(foreign)

# Inputs
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train.rds")
ytrain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_full_train.csv")
folds <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train_folds.csv")
output_location <- "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/neuro_full_train"
cv <- 12

# Sanity check 
identical(folds$Basename, ytrain$Basename)
identical(rownames(xtrain), ytrain$Basename)

# Load model
analysis <- function(xtrain, ytrain, cv, folds) {
 
  set.seed(1783) # Set seed
  seed <- 1783

  Output <- list() # Create storage for predictor weights results
  Output2 <- list() # Create storage for the features table summary 

  for (i in cols) { # Loop through the proteins by columns index which can be set below, before running the analysis function
    tryCatch({
      message("1. Assigning the protein variable.")
      x <- xtrain
      q <- ytrain[1] # Get just basenames for people in the y variable 
      p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
      name_p <- colnames(p) # Get the name of the protein for this iteration
      y <- cbind(q,p) # Bind Basename and protein data together into one set 
      index <- identical(y$Basename, rownames(x)) # check this has been done
      names(y)[2] <- "protein" # Assign a generic name to the protein variable
      y <- as.numeric(y$protein) # Create a numeric list for this variable to feed in as y to the model

      message(paste0("2. Computing predictor weights for ", name_p)) 
      fold <- folds$folds
      lasso.cv <- cv.glmnet(x, y, family="gaussian", alpha = 0.5, foldid = fold, nfolds=cv) # cross validation to get best lambda
      fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, standardize = F, lambda = lasso.cv$lambda.min) # model fit 
      coefs <- coef(fit) # Extract coeficients 
      coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
      coefs$Predictor <- name_p # Assign protein identifier
      names(coefs)[1] <- "Coefficient" # Tidy naming 
      coefs$CpG_Site <- rownames(coefs) # Create cpg column
      coefs <- coefs[c(3,1,2)] # order 
      coefs3 <- coefs[-1,] # Remove intercept (if there was only intercept, this will now have 0 rows)
      write.csv(coefs3, file = paste0(output_location, "/01_predictor_weights/predictor_weights_for_protein_", name_p, ".csv"), row.names = F)

      features <- nrow(coefs3) # Count number of features contributing for the protein - i.e. if its > 0 it has features 
      Statement <- matrix(nrow = (length(cols)), ncol = 2) %>% as.data.frame() # Create a vessel to store the results from the 

      if (features > 0) {
        rownames(coefs3) <- NULL # remove rows as they arent needed 
        results <- coefs3 # save weights as outputs list to bind together for all proteins
        message(paste0("3. Your weights for protein ", name_p, " have been calculated."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
        
      } else {
        results <- NULL
        message(paste0("3. Your model for protein ", name_p, " did not generate features."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
      }
      Output[[i]] <- results
      Output2[[i]] <- Statement
    }, error = function(e) cat("skipped"))
  }

  weight_total <- do.call(rbind, Output) # Bind predictor weights together
  write.csv(weight_total, paste0(output_location, "/weights_output_210121.csv"), row.names = F)

  statement_total <- do.call(rbind, Output2) # Bind summary table together 
  statement_total <- na.omit(statement_total)
  names(statement_total) <- c("Uniprot_ID", "Features")
  write.csv(statement_total, paste0(output_location, "/statement_output_210121.csv"), row.names = F)
  
  return(statement_total)
  message("Model running complete! The outputs will be in your specified location.")
}

message("6. Run the model.")

cols <- c(2:93) # Define cols to loop through here for the protein variables 
cv <- 12 # tell it which cv to run for 

start_time <- Sys.time()
results <- analysis(xtrain = xtrain, ytrain = ytrain, cv = cv, folds = folds)
end_time <- Sys.time()
time_difference <- end_time - start_time 




############################################################################################

### INFLAM TRAIN/TEST

# Required libraries
library(glmnet)
library(tidyverse)

# Inputs
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_train_only.rds")
xtest <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_test_only.rds")
xtest <- t(xtest)
ytrain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_train_only.csv")
ytest <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_test_only.csv")
folds <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_train_only_folds.csv")
output_location <- "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/inflam_train_test"
cv <- 10

# Sanity check 
identical(folds$Basename, ytrain$Basename)
identical(rownames(xtrain), ytrain$Basename)
identical(colnames(xtest), ytest$Basename)

# Load model
analysis <- function(xtrain, ytrain, xtest, ytest, cv, folds) {
 
  set.seed(1783) # Set seed
  seed <- 1783
  cv <- 10

  Output <- list() # Create storage for predictor weights results
  Output2 <- list() # Create storage for the features table summary 
  Output3 <- list() # Create storage for the predictor scores results

  for (i in cols) { # Loop through the proteins by columns index which can be set below, before running the analysis function
    tryCatch({
      message("1. Assigning the protein variable.")
      x <- xtrain
      q <- ytrain[1] # Get just basenames for people in the y variable 
      p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
      name_p <- colnames(p) # Get the name of the protein for this iteration
      y <- cbind(q,p) # Bind Basename and protein data together into one set 
      index <- identical(y$Basename, rownames(x)) # check this has been done
      names(y)[2] <- "protein" # Assign a generic name to the protein variable
      y <- as.numeric(y$protein) # Create a numeric list for this variable to feed in as y to the model

      message(paste0("2. Computing predictor weights for ", name_p)) 

      fold <- folds$folds
      lasso.cv <- cv.glmnet(x, y, family="gaussian", alpha = 0.5, foldid = fold, nfolds=cv) # cross validation to get best lambda
      fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, standardize = F, lambda = lasso.cv$lambda.min) # model fit 
      coefs <- coef(fit) # Extract coeficients 
      coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
      coefs$Predictor <- name_p # Assign protein identifier
      names(coefs)[1] <- "Coefficient" # Tidy naming 
      coefs$CpG_Site <- rownames(coefs) # Create cpg column
      coefs <- coefs[c(3,1,2)] # order 
      coefs3 <- coefs[-1,] # Remove intercept (if there was only intercept, this will now have 0 rows)

      features <- nrow(coefs3) # Count number of features contributing for the protein - i.e. if its > 0 it has features 
      Statement <- matrix(nrow = (length(cols)), ncol = 2) %>% as.data.frame() # Create a vessel to store the results

      if (features > 1) {

        message(paste0("3. Computing predictor scores for ", name_p)) 

        overlap <- which(rownames(xtest) %in% rownames(coefs3)) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
        meth <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
        match <- meth[match(rownames(coefs3), rownames(meth)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
        calc <- match * coefs3[,2] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
        sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
        export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
        names(export_sum)[1] <- "Scores"
        export_sum$Predictor <- name_p
        export_sum$Basenames <- rownames(export_sum)
        write.csv(export_sum, file = paste0(output_location, "/01_predictor_scores/predictor_scores_for_protein_", name_p, ".csv"), row.names = F)

        rownames(coefs3) <- NULL # remove rows as they arent needed 
        results <- coefs3 # save weights as outputs list to bind together for all proteins
        message(paste0("4. Your scores for protein ", name_p, " have been calculated."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
      } else {
        results <- NULL
        export_sum <- NULL
        message(paste0("4. Your model for protein ", name_p, " did not generate features."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
      }
      Output[[i]] <- results
      Output2[[i]] <- Statement
      Output3[[i]] <- export_sum
    }, error = function(e) cat("skipped"))
  }

  weight_total <- do.call(rbind, Output) # Bind predictor weights together
  write.csv(weight_total, paste0(output_location, "/weights_output_200121.csv"), row.names = F)

  statement_total <- do.call(rbind, Output2) # Bind summary table together 
  statement_total <- na.omit(statement_total)
  names(statement_total) <- c("Uniprot_ID", "features")
  write.csv(statement_total, paste0(output_location, "/statement_output_200121.csv"), row.names = F)

  score_total <- do.call(rbind, Output3) # Bind together scores 
  write.csv(score_total, paste0(output_location, "/score_output_200121.csv"), row.names = F)
  
  return(statement_total)
  message("Model running complete! The outputs will be in your specified location.")
}

message("6. Run the model.")

cols <- c(2:71) # Define cols to loop through here for the protein variables 2:71
cv <- 10 # tell it which cv to run for 

start_time <- Sys.time()
results <- analysis(xtrain = xtrain, ytrain = ytrain, xtest = xtest, ytest = ytest, cv = cv, folds = folds)
end_time <- Sys.time()
time_difference <- end_time - start_time 

####################################################################################################################

### ASSESSMENT: CORRELATE SCORES IN HOLDOUT SAMPLE

# The xtest set was used to calcualte scores above
# I will correlate the DNAm-generated scores from xtest with the ytest protein measurements
# I'll extract a results table with correlation, pval and CIs

table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/table_pQTLs.csv")

library(tidyverse)
ytest <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_test_only.csv")
scores <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/inflam_train_test/score_output_200121.csv")
summary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/inflam_train_test/statement_output_200121.csv")

# Look at how many generted features 
number <- summary %>% filter(features > 0) # 48 changed to 50 with pQTLs regressed

# Look at how many generated scores 
number <- scores$Predictor %>% unique() # 44 has remained the same at 44 with pQTLs regressed 

# Spread the scores so one protein per column
test <- spread(scores, Predictor, Scores)

# Check to make sure they have spread correctly 
ADA <- test[c(1:2)] # spread version
ADA2 <- scores %>% filter(Predictor == "ADA")
people <- ADA2$Basename
ADA_match <- ADA[match(people, ADA$Basename),] # match ADA to the original ADA from scores file which was unspread
identical(ADA_match$Basename, ADA2$Basename) # TRUE 
identical(ADA_match$ADA, ADA2$Scores) # TRUE - the scores match exactly 

inflam <- test

# Match the inflam file to the ytest file in terms of basenames
inflam <- inflam[match(ytest$Basename, inflam$Basenames),]

# Check this 
inflam$Basenames <- as.character(inflam$Basenames)
identical(ytest$Basename, inflam$Basenames)
names(inflam)[1] <- "Basename"

# Match order of columns between the 2 files for proteins with scores 
overlap <- which(colnames(ytest) %in% colnames(inflam))
ytest <- ytest[,overlap]
names <- colnames(ytest)
inflam <- inflam[match(names, colnames(inflam))]
identical(colnames(ytest), colnames(inflam)) # TRUE

index <- inflam$Basename

# Remove basenames
inflam <- inflam[-1]
ytest <- ytest[-1]

# Make results table 
results <- data.frame(Protein = colnames(inflam)[1:44], r = 1:44, p = 1:44, LC = 1:44, UC = 1:44)

# Correlate for each protein against the proxy 

for (i in 1:44){
name <- as.character(colnames(inflam)[i])
cor1 <- cor.test(ytest[,name], inflam[,name])
int <- cor1$conf.int[1:2]
p <- cor1$p.value 
r <- cor1$estimate 
results[i,"r"] <- r
results[i,"p"] <- p
results[i,4:5] <- int
}

# Order by r
results <- results[rev(order(results$r)),]

# Apply cutoff of r > 0.1 and p < 0.05 
result <- results %>% filter(r > 0.1)
result2 <- result %>% filter(p < 0.1)

# Write out the new results for comparison reference 
write.csv(result2, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/inflam_pass_new.csv", row.names = F)

############################################################################################

### NEURO TRAIN/TEST

# Required libraries
library(glmnet)
library(tidyverse)

# Inputs
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_train_only.rds")
xtest <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_test_only.rds")
xtest <- t(xtest)
ytrain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_train_only.csv")
ytest <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_test_only.csv")
folds <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_train_only_folds.csv")
output_location <- "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/neuro_train_test"
cv <- 10

# Sanity check 
identical(folds$Basename, ytrain$Basename)
identical(rownames(xtrain), ytrain$Basename)
identical(colnames(xtest), ytest$Basename)

# Load model
analysis <- function(xtrain, ytrain, xtest, ytest, cv, folds) {
 
  set.seed(1783) # Set seed
  seed <- 1783
  cv <- 10

  Output <- list() # Create storage for predictor weights results
  Output2 <- list() # Create storage for the features table summary 
  Output3 <- list() # Create storage for the predictor scores results

  for (i in cols) { # Loop through the proteins by columns index which can be set below, before running the analysis function
    tryCatch({
      message("1. Assigning the protein variable.")
      x <- xtrain
      q <- ytrain[1] # Get just basenames for people in the y variable 
      p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
      name_p <- colnames(p) # Get the name of the protein for this iteration
      y <- cbind(q,p) # Bind Basename and protein data together into one set 
      index <- identical(y$Basename, rownames(x)) # check this has been done
      names(y)[2] <- "protein" # Assign a generic name to the protein variable
      y <- as.numeric(y$protein) # Create a numeric list for this variable to feed in as y to the model

      message(paste0("2. Computing predictor weights for ", name_p)) 

      fold <- folds$folds
      lasso.cv <- cv.glmnet(x, y, family="gaussian", alpha = 0.5, foldid = fold, nfolds=cv) # cross validation to get best lambda
      fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, standardize = F, lambda = lasso.cv$lambda.min) # model fit 
      coefs <- coef(fit) # Extract coeficients 
      coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
      coefs$Predictor <- name_p # Assign protein identifier
      names(coefs)[1] <- "Coefficient" # Tidy naming 
      coefs$CpG_Site <- rownames(coefs) # Create cpg column
      coefs <- coefs[c(3,1,2)] # order 
      coefs3 <- coefs[-1,] # Remove intercept (if there was only intercept, this will now have 0 rows)

      features <- nrow(coefs3) # Count number of features contributing for the protein - i.e. if its > 0 it has features 
      Statement <- matrix(nrow = (length(cols)), ncol = 2) %>% as.data.frame() # Create a vessel to store the results

      if (features > 1) {

        message(paste0("3. Computing predictor scores for ", name_p)) 

        overlap <- which(rownames(xtest) %in% rownames(coefs3)) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
        meth <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
        match <- meth[match(rownames(coefs3), rownames(meth)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
        calc <- match * coefs3[,2] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
        sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
        export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
        names(export_sum)[1] <- "Scores"
        export_sum$Predictor <- name_p
        export_sum$Basenames <- rownames(export_sum)
        write.csv(export_sum, file = paste0(output_location, "/01_predictor_scores/predictor_scores_for_protein_", name_p, ".csv"), row.names = F)

        rownames(coefs3) <- NULL # remove rows as they arent needed 
        results <- coefs3 # save weights as outputs list to bind together for all proteins
        message(paste0("4. Your scores for protein ", name_p, " have been calculated."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
      } else {
        results <- NULL
        export_sum <- NULL
        message(paste0("4. Your model for protein ", name_p, " did not generate features."))
        Statement[1,1] <- name_p
        Statement[1,2] <- features
      }
      Output[[i]] <- results
      Output2[[i]] <- Statement
      Output3[[i]] <- export_sum
    }, error = function(e) cat("skipped"))
  }

  weight_total <- do.call(rbind, Output) # Bind predictor weights together
  write.csv(weight_total, paste0(output_location, "/weights_output_200121.csv"), row.names = F)

  statement_total <- do.call(rbind, Output2) # Bind summary table together 
  statement_total <- na.omit(statement_total)
  names(statement_total) <- c("Uniprot_ID", "features")
  write.csv(statement_total, paste0(output_location, "/statement_output_200121.csv"), row.names = F)

  score_total <- do.call(rbind, Output3) # Bind together scores 
  write.csv(score_total, paste0(output_location, "/score_output_200121.csv"), row.names = F)
  
  return(statement_total)
  message("Model running complete! The outputs will be in your specified location.")
}

message("6. Run the model.")

cols <- c(2:93) # Define cols to loop through here for the protein variables 
cv <- 10 # tell it which cv to run for 

start_time <- Sys.time()
results <- analysis(xtrain = xtrain, ytrain = ytrain, xtest = xtest, ytest = ytest, cv = cv, folds = folds)
end_time <- Sys.time()
time_difference <- end_time - start_time 

#############################################################################################################

### ASSESSMENT: CORRELATE SCORES IN HOLDOUT SAMPLE

# The xtest set was used to calcualte scores above
# I will correlate the DNAm-generated scores from xtest with the ytest protein measurements
# I'll extract a results table with correlation, pval and CIs

library(tidyverse)
ytest <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_test_only.csv")
scores <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/neuro_train_test/score_output_200121.csv")
summary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/neuro_train_test/statement_output_200121.csv")

# Look at how many generted features 
number <- summary %>% filter(features > 0) # 63 changed to 61 with features 

# Look at how many generated scores 
number <- scores$Predictor %>% unique() # 63 to 61 with scores 

# Spread the scores so one protein per column
test <- spread(scores, Predictor, Scores)

neuro <- test

# Match the file to the ytest file in terms of basenames
neuro <- neuro[match(ytest$Basename, neuro$Basenames),]

# Check this 
neuro$Basenames <- as.character(neuro$Basenames)
identical(ytest$Basename, neuro$Basenames)
names(neuro)[1] <- "Basename"

# Match order of columns between the 2 files for proteins with scores 
overlap <- which(colnames(ytest) %in% colnames(neuro))
ytest <- ytest[,overlap]
names <- colnames(ytest)
neuro <- neuro[match(names, colnames(neuro))]
identical(colnames(ytest), colnames(neuro)) # TRUE

index <- neuro$Basename

# Remove basenames
neuro <- neuro[-1]
ytest <- ytest[-1]

# Make results table 
results <- data.frame(Protein = colnames(neuro)[1:61], r = 1:61, p = 1:61, LC = 1:61, UC = 1:61)

# Correlate for each protein against the proxy 

for (i in 1:61){
name <- as.character(colnames(neuro)[i])
cor1 <- cor.test(ytest[,name], neuro[,name])
int <- cor1$conf.int[1:2]
p <- cor1$p.value 
r <- cor1$estimate 
results[i,"r"] <- r
results[i,"p"] <- p
results[i,4:5] <- int
}

# Order by r
results <- results[rev(order(results$r)),]

# Apply cutoff of r > 0.1 and p < 0.05 
result <- results %>% filter(r > 0.1)
result2 <- result %>% filter(p < 0.1)

# Write out the new results for comparison reference 
write.csv(result2, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Comparison_model_2_output/neuro_pass_new.csv", row.names = F)
