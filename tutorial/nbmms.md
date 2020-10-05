# NBZIMM - NBMM (Negative Binomial Mixed Model)

## Introduction

The complex microbiome is inherently dynamic. The metagenomics sequencing data provide valuable resources for investigating the dynamic changes of microbial abundance over time and the associations between the microbiome and host environmental/clinical factors. The well-known properties of microbiome measurements include varied total sequence reads across samples, over-dispersion and zero-inflation. Additionally, microbiome studies usually collect samples longitudinally, which insert correlation among the samples and thus further complicate the analysis and interpretation of microbiome count data. In this tutorial, we implement our proposed Negative Binomial mixed models (NBMMs) for detecting the association between the microbiome and host environmental/clinical factors for longitudinal microbiome data.

## Usage
```r
glmm.nb(fixed, random, data, subset, correlation, weights, control, niter = 30, epsilon = 1e-05, verbose = TRUE)
```

## Arguments

- **fixed, random, data, subset, correlation, weights, control**: These arguments are the same as in the function lme in the package nlme.
- **niter**: maximum number of iterations.
- **epsilon**: positive convergence tolerance.
- **verbose**: logical. If TRUE, print out number of iterations and computational time.

## Real Data Analysis Examples
The following R code is used for real data analysis in a manuscript and the citation will added later.

You can find the datasets: <https://github.com/abbyyan3/NBZIMM-tutorial/tree/master/NBMM-longitudinal-temporal-data/>

```r
library(phyloseq)
library(NBZIMM)
library(nlme)

setwd("directory")

clinical <- read.csv("temporal spatial clinical merged.csv")
otu <- read.csv("temporal spatial otu.csv", check.names = F)
taxonomy <- read.csv("temporal spatial taxonomy.csv")
clusters <- read.csv("clust.csv")

clinical <- merge(clinical, clusters, by.x = "X", by.y = "X")
########### otu data #################
dim(otu)
colnames(otu)
rownames((otu))

############### clinical data ############################
dim(clinical)
colnames(clinical)

sites = clinical$BodySite  # four different body sites
table(sites)

# Samples collected before delivery
Preg <- clinical$Preg
table(sites, Preg)
PrePreg <- clinical$PrePreg
table(Preg, PrePreg)

# Samples collected after delivery
Post <- clinical$PostPreg

#### keep only samples collected during pregnancy
##### Keep only samples in vaginal swab ####
### To compare with CST analysis in the paper
clinical <- clinical[clinical$Outcome != "Marginal", ] ### Same as in the paper
clinical <- clinical[clinical$V1 == "4", ] ### Same as in the paper
clinical <- clinical[clinical$Preg == "TRUE", ]
clinical <- clinical[clinical$BodySite == "Vaginal_Swab", ] 
otu <- otu[otu$Sample_ID %in% clinical$SampleID, ]

# Remove taxa with proportions of zeros greater than 0.25
yy = as.matrix(otu[,-1])
yy = ifelse(is.na(yy), 0, yy)
zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.25, ]) ## same as the proportion of non zero used in the paper for CST = 4
yy = yy[, rownames(zero.p)]
rownames(yy) = otu$Sample_ID

sampleID = clinical$SampleID # different samples for one subject
outcome = clinical$Outcome # preterm vs term groups

# Variables for fitting the models
SubjectID = as.factor(clinical$SubjectID) # subject id for random effects 
term_label = ifelse(outcome == "Term", 0, 1) # group indicator
weeks = as.numeric(log(clinical$GDColl)) # Days for samples collected
N = clinical$NumReads  # total reads for library sizes


##########################
## To fit a NBMM model, we need to use glmm.nb funtion to analyze the taxa in raw counts.
### In the following code, we compared four different models using the function glmm.nb to separately analyze each taxon for all taxa in the real data.

# choose the first taxon to demonstrate the examples:
    y = as.numeric(yy[, 1])
 
 # The first model compared each taxon between two groups and an offset term for library size is needed in fixed effects. The random effect includes a random intercept.
f1 = glmm.nb(y ~ term_label + offset(log(N)), random = ~ 1|SubjectID)
# summarize the results
summary(f1)

 # The second model compared each taxon between two groups and an offset term for library size is needed in fixed effects. The random effects include a random intercept and a random slope with weeks.
f2 = glmm.nb(y ~ term_label + offset(log(N)), random = list(SubjectID = pdDiag(~weeks)))
summary(f2)

 # The third model compared each taxon between two groups, while controlling covariate weeks and interaction between group and weeks. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept.
f3 = glmm.nb(y ~ term_label*weeks + offset(log(N)), random = ~ 1|SubjectID)
summary(f3)

 # The fourth model compared each taxon between two groups, while controlling covariate weeks and interaction between group and weeks. And an offset term for library size is needed in fixed effects. The random effects include a random intercept and a random slope with weeks.
f4 = glmm.nb(y ~ term_label*weeks + offset(log(N)), random = list(SubjectID = pdDiag(~weeks)))
summary(f4)

For all taxa, we can analyze them using the fuction mms for all four above models:

# The first model:
f = mms(y = yy, fixed = ~ term_label + offset(log(N)), 
        random = ~ 1 | SubjectID, data = data,
        min.p = 0.2, method = "nb")

# The second model:
f = mms(y = yy, fixed = ~ term_label + offset(log(N)), 
        random =list(SubjectID = pdDiag(~weeks)), data = data,
        min.p = 0.2, method = "nb")
        
# The third model:
f = mms(y = yy, fixed = ~ term_label*weeks + offset(log(N)), 
        random = ~ 1 | SubjectID, data = data,
        min.p = 0.2, method = "nb")
        
# The fourth model:
f = mms(y = yy, fixed = ~ term_label*weeks + offset(log(N)), 
        random = list(SubjectID = pdDiag(~weeks)), data = data,
        min.p = 0.2, method = "nb")
       
```

