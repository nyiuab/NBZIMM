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


for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    f17 = glmm.zinb(y ~ group + offset(log(N))|1, random = ~ 1|subject)
    out29[j, ] = summary(f17)$tTable[2, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f18 = glmm.zinb(y ~ group*days + offset(log(N))|1, random = list(subject = pdDiag(~days)))
    out30[j, ] = summary(f18)$tTable[2, 5]
    out31[j, ] = summary(f18)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f19 = glmm.zinb(y ~ group*days + offset(log(N))|1, random = ~ 1|subject)
    out32[j, ] = summary(f19)$tTable[2, 5]
    out33[j, ] = summary(f19)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f20 = glmm.zinb(y ~ group*days  + offset(log(N))|1, random = ~ 1|subject, correlation = corAR1()) #random = list(subject = pdDiag(~days)))
    out34[j, ] = summary(f20)$tTable[2, 5]
    out35[j, ] = summary(f20)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f21 = glmm.zinb(y ~ group + offset(log(N)), random = ~ 1|subject)
    out36[j, ] = summary(f21)$tTable[2, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f22 = glmm.zinb(y ~ group*days + offset(log(N)), random = list(subject = pdDiag(~days)))
    out37[j, ] = summary(f22)$tTable[2, 5]
    out38[j, ] = summary(f22)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f23 = glmm.zinb(y ~ group*days + offset(log(N)), random = ~ 1|subject)
    out39[j, ] = summary(f23)$tTable[2, 5]
    out40[j, ] = summary(f23)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    
    f24 = glmm.zinb(y ~ group*days  + offset(log(N)),random = ~ 1|subject, correlation = corAR1()) # random = list(subject = pdDiag(~days)))
    out41[j, ] = summary(f24)$tTable[2, 5]
    out42[j, ] = summary(f24)$tTable[4, 5]      
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
```
