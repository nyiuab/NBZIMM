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
library(NBZIMM)

data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)

N = sam[, "Total.Read.Counts"]        
Days = sam$GA_Days; Days = scale(Days)
Age = sam$Age; Age = scale(Age)
Race = sam$Race
preg = sam$pregnant; table(preg)
subject = sam[, "Subect_ID"]; table(subject)

y = otu[, 1]

f1 = glmm.nb(y ~ Days + Age + Race + preg + offset(log(N)), random = ~ 1 | subject)
summary(f1)
fixed(f1)

For all taxa, we can analyze them using the fuction mms for all four above models:

# The first model:
f2 = mms(y = Romero$OTU, fixed = ~  Days + Age + Race + preg + offset(log(N)), 
        random = ~ 1 | subject, data = Romero$SampleData,
        min.p = 0.2, method = "nb")


       
```

