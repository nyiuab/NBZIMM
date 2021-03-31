# NBZIMM - ZINBMMs (Zero-inflated Negative Binomial Mixed Model)

## Introduction

This function sets up and fits zero-inflated negative binomial mixed models for analyzing zero-inflated count responses with multilevel data structures (for example, clustered data and longitudinal studies).

## Usage
```r
glmm.zinb(fixed, random, data, correlation, zi.random = FALSE, 
        niter = 30, epsilon = 1e-05, verbose = TRUE, ...) 
```
## Arguments

- **fixed**: symbolic description of the fixed-effects part of the model, see details.
- **random, data, correlation**: These arguments are the same as in the function lme in the package nlme.
- **zi.random**: logical. If TRUE, include the random effect specified by random in the zero-inflation part.
- **niter**: maximum number of iterations.
- **epsilon**: positive convergence tolerance.
- **verbose**: logical. If TRUE, print out number of iterations and computational time.
...**: further arguments for lme.

## Real Data Analysis Examples:
The following R code is used for real data analysis in a manuscript and the citation will added later.
You can find the datasets: <https://github.com/abbyyan3/NBZIMM-tutorial/tree/master/ZINBMMs-dataset>

The Romero's preterm DATASET
```r
rm(list = ls())

library(NBZIMM)
library(nlme)

load("longitudinal pregnant data term vs preterm.RData")
########### otu data #################
dim(otu)
colnames(otu)

############### clinical data ############################
clinical = data.frame(clinical)
dim(clinical)
colnames(clinical)

group = ifelse(clinical$Pregnancy.outcome == "TERM", 1, 0)

# Library size as offset
N = clinical[, "Total.Read.Countsb"]  # total reads

# Variables for fitting the models
subject = clinical[, "subject.ID"]
clinical$group = group
## pregnant women vs non-pregnant women measurement time in weeks
days = clinical[, "calculated.GA..weeks."]

# Remove taxa with proportions of zeros between 0.2-0.8.
pheno = round(otu[, -c(1,2)]/100*N)
yy = as.matrix(pheno)  
yy = ifelse(is.na(yy), 0, yy)
zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.2 & zero.p$zero.p<0.8, ])
yy = yy[, rownames(zero.p)]

##########################
## To fit a ZINBMM model, we need to use glmm.zinb funtion to analyze the taxa count data.
### In the following code, we compared five different models using the function glmm.zinb to separately analyze each taxon for all taxa in the real data.

# choose the first taxon to demonstrate the examples:
    y = as.numeric(yy[, 1])

data = cbind.data.frame(y, subject, group, N)

# The first model compared each taxon between two groups and an offset term for library size is needed in fixed effects. The random effect includes a random intercept. No fixed or random effects is considered for zero part. 
    f1 = glmm.zinb(y ~ group + offset(log(N)), data = data, 
               random = ~ 1|subject, zi_fixed = ~1, zi_random =NULL)
    summary(f1)

# The second model compared each taxon between two groups, while controlling days, and group by days interaction. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept and a random slope of days. No fixed or random effects is considered for zero part. 
    f2 = glmm.zinb(y ~ group*days + offset(log(N)), data = data, 
               random = list(subject = pdDiag(~days)), zi_fixed = ~1, zi_random =NULL)
    summary(f2)
    
# The third model compared each taxon between two groups, while controlling days, and group by days interaction. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept. No fixed or random effects is considered for zero part.
    f3 =  glmm.zinb(y ~ group*days + offset(log(N)), data = data, 
               random = ~ 1|subject, zi_fixed = ~1, zi_random =NULL)
    summary(f3)

# The fourth model compared each taxon between two groups, while controlling days, and group by days interaction. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept. A correlation structure of AR1 is considered in this model. No fixed or random effects is considered for zero part.
    f4 = glmm.zinb(y ~ group*days  + offset(log(N)), data = data, 
               random = ~ 1|subject, zi_fixed = ~1, zi_random =NULL, correlation = corAR1()) 
    summary(f4)
    
# The fifth model compared each taxon between two groups and an offset term for library size is needed in fixed effects. The random effect includes a random intercept. Group is considered for fixed effects in zero part. 
    f5 = glmm.zinb(y ~ group + offset(log(N)), zi_fixed = ~group, zi_random =NULL,, random = ~ 1|subject)
    summary(f5)
       
For all taxa, we can analyze them using the fuction mms for all the above models. For demonstration, we use mms for the first model:
# The first model:
f = mms(y = yy, fixed = group + offset(log(N)), 
        random = ~ 1 | subject, data = data,zi_fixed = ~1, zi_random =NULL,
        min.p = 0.2, method = "zinb")
        
```
