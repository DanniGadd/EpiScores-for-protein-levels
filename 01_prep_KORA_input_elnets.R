Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# INPUT PREP KORA TRAINING #################################################
############################################################################################
############################################################################################

# This script preps KORA input files for elnet models 

# libraries
library(tidyverse)
library(glmnet)
library(imputeTS)
library(bestNormalize)
library(foreign)

###################################
# load KORA data

#load proteins 
load("PROTEINSraw.RData")

#read the 793 common proteins between KORA and STRADL
seqIDs_KORA_STRADL_793 <- read.csv("seqIDs_KORA_STRADL_793.csv", sep="")

#match the 793 common proteins in KORA (using the shortened SeqIDs)
PROTEINSraw=PROTEINSraw[match(as.matrix(seqIDs_KORA_STRADL_793),sapply(strsplit(rownames(PROTEINSraw),"_") , "[[", 1)),]

#load methylation M-values
load("mvalue.RData")

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}
betas <- m2beta(mvalue) # convert m-value to beta-value 


#read the protein annotation
ProteinSeqID_positions <- read.delim("ProteinSeqID_positions.txt", dec=",")


#load covariates
load("covariates.RData")
age=covariates$age
sex=covariates$sex
pc1=covariates$PC1
pc2=covariates$PC2
pc3=covariates$PC3
pc4=covariates$PC4
pc5=covariates$PC5
pc6=covariates$PC6
pc7=covariates$PC7
pc8=covariates$PC8
pc9=covariates$PC9
pc10=covariates$PC10
pc11=covariates$PC11
pc12=covariates$PC12
pc13=covariates$PC13
pc14=covariates$PC14
pc15=covariates$PC15
pc16=covariates$PC16
pc17=covariates$PC17
pc18=covariates$PC18
pc19=covariates$PC19
pc20=covariates$PC20


#load the imputed KORA SNP data
load("SNP_annotation_imputed_renamed.RData")
load("KORAsnps_imputed_renamed.RData")
rsid_snp=as.matrix(SNP_annotation$rsID)

id_table <- read.delim("id_table.txt")

#match the SNP ids to the CpG ids
rownames(KORAsnps)=id_table$cpg_id[match(rownames(KORAsnps),id_table$snp_id)]

#ensure the renaming was done correctly
identical(rownames(KORAsnps),colnames(PROTEINSraw))

fvalid = 0.1 # fraction of data to be set aside for validation

# define the proteins to be tested here
range = c(1:793)

# fix random seed to make things reproducible
set.seed(4711)

# the methylation data (x) must be an array NSAMPLES x NCPG
# rownames are sample ids 
# colnames are cpg ids
#x = t(mvalue)
x = t(betas)

## preprocess methylation data 
## Replace NA with imputed values with CpGs as columns 
#number of NA before imputation
length(which(is.na(x)==TRUE))

#number of NA before impuation
x <- na_mean(x)

#number of NA after imputation
length(which(is.na(x)==TRUE))

# Scale the methylation dataset with CpGs as columns
zwi = rownames(x)
x <- apply(x, 2, scale)
rownames(x) <- zwi

# set a fraction of the data aside as validation set
ixvalid = sample(seq(dim(x)[1]),fvalid*dim(x)[1])
xvalid = x[ixvalid,]
#x = x[-ixvalid,]

# log transform protein data
logPro=log(PROTEINSraw)
d<-t(logPro)

#residualize
d=as.data.frame(d)
for(i in ncol(d)) { # Residualise proteins by desired phenotypes and scale 
  d[,i] <- scale(resid(lm(d[,i] ~ age + sex  + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + 
  pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + pc16 + pc17 + pc18 + pc19 +pc20, 
  data=d, na.action="na.exclude")))
}

###Correct proteins for SNPs#############
#read the protein annotation
ProteinSeqID_positions <- read.delim("ProteinSeqID_positions.txt", dec=",")
pro_Target=as.matrix(ProteinSeqID_positions$SeqId)
pro_Target_chr=as.matrix(ProteinSeqID_positions$chromosome)
pro_Target_pos_start=as.matrix(ProteinSeqID_positions$start_position)
pro_Target_pos_end=as.matrix(ProteinSeqID_positions$end_position)
pro_Target_EntrezGeneSymbol=as.matrix(ProteinSeqID_positions$EntrezGeneSymbol)
pro_Target_Target=as.matrix(ProteinSeqID_positions$Target)
pro_Target_TargetFullName=as.matrix(ProteinSeqID_positions$TargetFullName)

library("readxl")
pGWAS_allhits <- read_excel("pGWAS_allhits.xlsx",  sheet = "Supplemental Data 5 ")

ind_matching=match(pGWAS_allhits$SNP,colnames(KORAsnps))
pGWAS_allhits=pGWAS_allhits[-which(is.na(ind_matching)==TRUE),]
ind_matching=ind_matching[-which(is.na(ind_matching)==TRUE)]
KORAsnps_from_PGWAS= KORAsnps[,ind_matching]
rsids_from_PGWAS=rsid_snp[ind_matching]  
  
library(stringr)
# It will automatically create a log file.
file.remove("Protein_SNP_assoc_LogFile.txt")
file.create("Protein_SNP_assoc_LogFile.txt")

for(i in 1:793)
{
  print(i)
  cat("-----------------------------\n",file="Protein_SNP_assoc_LogFile.txt",append=TRUE)
  cat("Protein :", i,'\n',file="Protein_SNP_assoc_LogFile.txt",append=TRUE)
  
  #get the list of SNPs for the protein
  snps_for_that_protein=pGWAS_allhits$SNP[which(pGWAS_allhits$SeqId ==  rownames(PROTEINSraw)[i])]
  #get the genotyping data for these SNPs
  genotyping_for_that_protein=KORAsnps[,match(snps_for_that_protein,rsid_snp)]
  genotyping_for_that_protein=as.matrix(genotyping_for_that_protein)
  
  if(length(genotyping_for_that_protein)>0)
{
        correct=TRUE
        while(correct)
        {
          tmp<-apply(genotyping_for_that_protein, 2, function(genotyping_for_that_protein){summary(lm(d[,i] ~ genotyping_for_that_protein ))$r.squared})
          cat("variance in snp ", colnames(genotyping_for_that_protein)[which(tmp==max(tmp))[1]],"is:",max(tmp),'\n',file="Protein_SNP_assoc_LogFile.txt",append=TRUE)
          #if that variance explained by that top SNP is less than 1% do nothing and exit 
          if(max(tmp)<0.01)
          {
            correct=FALSE
            cat("did not correct for any more SNPs\n",file="Protein_SNP_assoc_LogFile.txt",append=TRUE)
            break;
          }
          #else then correct for that SNP 
          else
          {
            m=lm(d[,i]~genotyping_for_that_protein[,which(tmp==max(tmp))[1]],na.action=na.exclude)
            #z-score the residuals
            r=scale(resid(m))
            #av=mean(r,na.rm=TRUE)
            #std=sd(r,na.rm=TRUE)
            #r=(r-av)/std
            d[,i]=r
            cat("corrected for SNP and z-scored\n",file="Protein_SNP_assoc_LogFile.txt",append=TRUE)
          }
        }
        
  
}
}

######################################################

## Rank-Inverse Based Normaliation
for(i in ncol(d)) {
  d[,i] <- orderNorm(d[,i])$x.t
}

## Scale the proteins with proteins as columns
zwi = rownames(d)
d <- apply(d, 2, scale)
rownames(d) <- zwi

#y = d[-ixvalid,]
y=d
yvalid=d[ixvalid,]


# subset to 397,630 CpG sites which are common to STRADL and GS cohorts,
d = read.csv("CpGs_sites_for_subsetting.csv")
shared_cpg = d[[1]]
shared_cpg=shared_cpg[-which(substr(shared_cpg,1,2)!='cg')]
diff_cpg = setdiff(shared_cpg,colnames(x))
shared_cpg = intersect(shared_cpg,colnames(x))
cat("there are", length(shared_cpg), "shared cpg sites\n")
cat("there are", length(diff_cpg), "cpgs that are not in KORA:\n")
cat("there are", length(which(substr(shared_cpg,1,2)!='cg')), "non-cg CpGs in the shared_cpg list\n")
head(diff_cpg) %>% print()
tail(diff_cpg) %>% print()

zwi = left_join(
  data.frame(id=shared_cpg),
  data.frame(id=colnames(x), ix = seq(length(colnames(x))))
)
ix = zwi$ix

cat("subsetting x\n")
x = x[,ix]
xvalid = xvalid[,ix]

save(x,y,xvalid,yvalid,ProteinSeqID_positions,file="allData.RData")




