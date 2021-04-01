
# NBZIMM - ZIGMMs (Zero-inflated Gaussian Mixed Model)

## Introduction

This function sets up and fits zero-inflated gaussian mixed models for analyzing zero-inflated continuous or count responses with multilevel data structures (for example, clustered data and longitudinal studies).

## Usage
```r
lme.zig(fixed, random, data, correlation, 
        zi_fixed = ~1, zi_random = NULL,
        niter = 30, epsilon = 1e-05, verbose = TRUE, ...)  
```
## Arguments

- **fixed**: a formula for the fixed-effects part of the Gaussian model, including the continuous normal outcome. This argument is the same as in the function lme in the package nlme.
- **random, data, subset, correlation**: These arguments are the same as in the function lme in the package nlme.
- **zi_fixed, zi_random**: formulas for the fixed and random effects of the zero inflated part. only contain the right-hand side part
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

data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)

N = sam[, "Total.Read.Counts"]  # total reads
preg = sam$pregnant; table(preg)
subject = sam[, "Subect_ID"]; table(subject)

non = nonzero(y = otu, total = N, plot = F)
nonzero.p = non[[1]]

N = sam[, "Total.Read.Counts"]        
Days = sam$GA_Days; Days = scale(Days)
Age = sam$Age; Age = scale(Age)
Race = sam$Race
preg = sam$pregnant; table(preg)
subject = sam[, "Subect_ID"]; table(subject)

y = otu[, 1]

y0 = log(y+1)
data = data.frame(y0=y0, Days=Days, Age=Age, Race=Race, preg=preg, N=N, subject=subject)
f1 = lme.zig(fixed = y0 ~ Days + Age + Race + preg + offset(log(N)), 
            random = ~ 1 | subject, data = data) 
summary(f1)
fixed(f1)
summary(f1$fit.zero)


f2 = mms(y = log(Romero$OTU+1), fixed = ~ GA_Days + Age + Race + pregnant + 
           offset(log(Total.Read.Counts)), data = Romero$SampleData,
         random = ~ 1 | subject, min.p = 0.2, method = "zig", 
         zi_fixed = ~1, zi_random = NULL)

```
