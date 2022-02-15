Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# Elastic net models for KORA ##############################################
############################################################################################
############################################################################################

# Here, we train elastic net models for 793 KORA proteins 

### INFLAM TRAIN FULL

# libraries
library(tidyverse)
library(glmnet)
library(imputeTS)
library(bestNormalize)
library(foreign)


load("allData.RData")
# open a tsv file

nproc=1

basename = "fit_KORA1"
outfile = paste(basename, nproc, "tsv", sep = ".")
cat ("output to", outfile, "\n")

Output = list()
k = 1
l = 1
cat(file = outfile, append = FALSE, "SeqId", "r", "p-value", "r valid", "p-value valid", "fit$df", "fit$dev.ratio", "fit$lambda", "label", sep = "\t", "\n")
#this next loop needs to be parallelized, otherwise will take ages to run
for (i in 1:793) {
  
  cat("working on protein", l, "of", length(range), "\n")
  
  ## Residualise proteins and scale - use covars
  odata = y[,i]
  SeqId = colnames(y)[i]
  label = ProteinSeqID_positions$EntrezGeneSymbol[match(colnames(y)[i],ProteinSeqID_positions$SeqId)]
  cat(SeqId, label, "\n")
  name_p = paste(SeqId, ":", label)
  
  
  
  start_time <- Sys.time()
  lasso.cv <- cv.glmnet(x, odata, family="gaussian", alpha = 0.5, nfolds=10) # cross validation to get best lambda (this would normally be set to 10)
  fit <- glmnet(x, odata, family = "gaussian", alpha = 0.5, standardize = F, lambda = lasso.cv$lambda.min) # model fit 
  save(fit,file=paste("fit",i,".RData",sep=""))
  end_time <- Sys.time()
  print( end_time - start_time)
  
  coefs <- coef(fit) # Extract coeficients 
  # coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  # coefs$Predictor <- name_p # Assign protein identifier
  # names(coefs)[1] <- "Coefficient" # Tidy naming 
  # coefs$CpG_Site <- rownames(coefs) # Create cpg column
  # coefs <- coefs[c(3,1,2)] # order 
  # coefs3 <- coefs[-1,] # Remove intercept (if there was only intercept, this will now have 0 rows)
  
  coefs3 = data.frame(
    CpG_site = coefs@Dimnames[[1]][coefs@i+1],
    Coefficient = coefs@x,
    Predictor = name_p
  )
  
  features <- nrow(coefs3) # Count number of features contributing for the protein - i.e. if its > 0 it has features 
  if (features > 0) {
    rownames(coefs3) <- NULL # remove rows as they arent needed 
    Output[[k]] <- coefs3 # save weights as outputs list to bind together for all proteins
    
    # test the prediction
    ypred = predict(fit, newx = x)
    c = cor.test(odata,ypred)
    
    yvalidpred = predict(fit, newx = xvalid)
    cvalid = cor.test(yvalid[,i],yvalidpred)
    
    cat(file = outfile, append = TRUE, SeqId, c$estimate, c$p.value, cvalid$estimate, cvalid$p.value, fit$df, fit$dev.ratio, fit$lambda, label, sep = "\t", "\n")
    cat(SeqId, c$estimate, cvalid$estimate, fit$df, label, sep = "\t", "\n")
    
    k = k+1
  } else {
    cat(file = outfile, append = TRUE, SeqId, rep(".",5), label, sep = "\t", "\n")
  }
  message(paste0(features, " weights for protein ", name_p, " have been calculated."))
  l = l+1
}

# save Output
weightfile = paste("weights", nproc, "csv", sep=".")
weight_total <- do.call(rbind, Output) # Bind predictor weights together
write.csv(weight_total, weightfile, row.names = F)
