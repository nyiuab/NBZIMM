
# NBZIMM - ZIGMMs (Zero-inflated Gaussian Mixed Model)

## Introduction

This function sets up and fits zero-inflated gaussian mixed models for analyzing zero-inflated continuous or count responses with multilevel data structures (for example, clustered data and longitudinal studies).

## Usage
```r
lme.zig(fixed, random, data, 
        zi.random = FALSE, correlation, niter = 30, epsilon = 1e-05, 
        verbose = TRUE, ...) 
```
## Arguments

- **fixed**: symbolic description of the fixed-effects part of the model, see details.
- **random, data, subset, correlation**: These arguments are the same as in the function lme in the package nlme.
- **zi.random**: logical. If TRUE, include the random effect specified by random in the zero-inflation part.
- **niter**: maximum number of iterations.
- **epsilon**: positive convergence tolerance.
- **verbose**: logical. If TRUE, print out number of iterations and computational time.
...**: further arguments for lme.

## Real Data Analysis
The following R code is used for real data analysis in a manuscript and the citation will added later.

The Romero's DATASET
```r
library(BhGLM)
library(NBZIMM)
library(nlme)
rm(list = ls())

data(Romero) ## data included in R package NBZIMM
names(Romero)
otu = Romero$OTU; dim(otu)
clinical = Romero$SampleData; dim(clinical)
colnames(clinical)
clinical = clinical[complete.cases(clinical[,3]), ]

group = c(rep(0, 139), rep(1, 758)) ### 0 for pregnant women vs 1 for non-pregnant women

# Library size as offset
N = clinical[, "Total.Read.Counts"]  # total reads

# Variables for fitting the models
subject = clinical[, "Subect_ID"]
age = clinical[, "Age"]
race = clinical[, "Race"]
clinical$group = group
## pregnant women vs non-pregnant women measurement time in weeks
weeks = ifelse(clinical$group == 0, clinical$GA_Days, clinical$GA_Days/7)

# Remove taxa with proportions of zeros between 0.1-0.3.
pheno = otu
yy = as.matrix(pheno)  
yy = ifelse(is.na(yy), 0, yy)
zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.1 & zero.p$zero.p<0.3, ])
yy = yy[, rownames(zero.p)]

##########################
## To fit a ZIGMM model, we need to use lme.zig funtion to analyze the transformed taxa data.
### In the following code, we compared four different models using the function lme.zig to separately analyze each taxon for all taxa in the real data.

# choose the first taxon to demonstrate the examples:
    y = as.numeric(yy[, 1])
    y1 = log2(y + 1) # use logrithmic transformation for count data
    
# The first model compared each taxon between two groups, while controlling weeks, age and race. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept. No fixed or random effects is considered for zero part. 
    f1 = lme.zig(y1 ~ group + weeks  + age + race + offset(log(N))|1, family = "zig", random = ~ 1|subject)
    summary(f1)

# The second model compared each taxon between two groups, while controlling weeks, group by weeks interaction, age and race. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept. No fixed or random effects is considered for zero part. 
    f2 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N))|1, family = "zig", random = ~ 1|subject)
    summary(f2)
    
# The third model compared each taxon between two groups, while controlling weeks, group by weeks interaction, age and race. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept. A correlation structure of AR1 is considered in this model. No fixed or random effects is considered for zero part.
    f3 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N))|1, family = "zig", random = ~ 1|subject, correlation = corAR1())
    summary(f3)

# The fourth model compared each taxon between two groups, while controlling weeks, group by weeks interaction, age and race. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept and a random slope of weeks. No fixed or random effects is considered for zero part.
    f4 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N))|1, family = "zig", random = list(subject = pdDiag(~weeks)))
    summary(f4)
    
# The fifth model compared each taxon between two groups, while controlling weeks, age and race. And an offset term for library size is needed in fixed effects. The random effect includes a random intercept. Group is considered for fixed effects in zero part. 
    f5 = lme.zig(y1 ~ group + weeks  + age + race + offset(log(N))|group, family = "zig", random = ~ 1|subject)
    summary(f5)
       
For all taxa, we can analyze them using the fuction mms for all the above models. For demonstration, we use mms for the first model:
# The first model:
f = mms(y = yy, fixed = group + weeks  + age + race + offset(log(N))|1, 
        random = ~ 1 | subject, data = data,
        min.p = 0.2, method = "zig")
```
