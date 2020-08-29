---
title: NBZIMM ZIGMMs Longitudinal
nav: true
---

# NBZIMM - ZIGMMs (Zero-inflated Gaussian Mixed Model)

## Introduction

This function sets up and fits zero-inflated gaussian mixed models for analyzing zero-inflated continuous or count responses with multilevel data structures (for example, clustered data and longitudinal studies).

## Installation
You can install our NBZIMM package by downloading NBZIMM_1.0.zip.
```r
install.packages("NBZIMM")
library(NBZIMM)
```

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

## Simulation Studies
The following R code is used for simulation studies (Setting 1) in a manuscript and the citation will added later.

```r
library(NBZIMM)
library(nlme)

### set number of individuals
n_2 <- 50
n_3 <- 100 # 
n_4 <- 150

n.ind <- n_2

n.t = 5              # number of individuals at each family or repeated measures
n = n.ind * n.t      # total number of individuals

a = rbind(c(0, 0), c(0.2, 0.3), c(0.3, 0.4))
b1 = a[1, ]
b2 = a[2, ]
b3 = a[3, ]


a1 = b1                    

theta0 = c(0.1, 5)      
mu0 = c(0.1, 3.5)        
sigma0 = c(0.5, 1)

p1 = c(0.0, 0.2)           
p2 = c(0.2, 0.4)           
p3 = c(0.4, 0.6)
p = p1

n.sims = 100
out1 = out2 = out3 = out4 = out5 = out6 = out7 = out8 = out9 = matrix(NA, n.sims, 3)

start.time <- Sys.time()
for (i in 1:n.sims){
  tryCatch({
    #set.seed(12345+i)    
    d = c(rep(0, n/2), rep(1, n/2))     # a binary variable: disease vs. healthy
    ind =  Time = rep(0, length(d))
    nd = 0
    for (j in 1:n.ind)
      for (k in 1:n.t) {
        nd = nd + 1
        ind[nd] = j   
        Time[nd] = log(k)
      } 
    
    sigma = runif(1, sigma0[1], sigma0[2]); sigma
    theta = runif(1, theta0[1], theta0[2]); theta
    p.zero = runif(1, p[1], p[2]); p.zero
    
    beta1 = runif(1, a1[1], a1[2])     # effect of d
    
    x = scale(d)
    
    ys = sim (n.ind = n.ind, n.measure = n.t, x.d = x, coef.d = beta1, tau.d = sigma, tau.z = 0,
              theta = theta, p.zero = p.zero)
    ind = ys$ind.ID
    x = ys$x.d
    y0 = ys$y
    off = log(ys$T)
    N = ys$T
    
    y2 = asin((y0/N)^0.5)
    f1 = lme(y2 ~ x, random = ~ 1 | ind) 
    out1[i, ] = summary(f1)$tTable[2, c(1,2,5)]
    
    y1 = log2((y0 + 1))
    f2 = lme.zig(y1 ~ x + offset(off)|1 , zi.random = F, random = ~ 1 | ind, verbose = F) 
    out4[i, ] = summary(f2)$tTable[2, c(1,2,5)]
    
    f4 = glmm.nb(y0 ~ x + offset(off), random = ~ 1 | ind, verbose = F) 
    out6[i, ] = summary(f4)$tTable[2, c(1,2,5)]    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
}
stop.time <- Sys.time()
round(difftime(stop.time, start.time, units = "min"), 3)

```

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

N = clinical[, "Total.Read.Counts"]  # total reads
mean(N); sd(N)
mean(log(N)); sd(log(N))

subject = clinical[, "Subect_ID"]
sample_id = clinical[, "Sample_ID"]
age = clinical[, "Age"]
race = clinical[, "Race"]
nugent_score = clinical[, "Nugent.Score"]
summary(nugent_score)
CST = clinical[, "CST"]
summary(CST)

clinical$group = group

## pregnant women vs non-pregnant women measurement time in weeks
weeks = ifelse(clinical$group == 0, clinical$GA_Days, clinical$GA_Days/7)

pheno = otu

pheno = pheno[rownames(clinical), ]
dim(pheno)
colnames(pheno)

yy = as.matrix(pheno)  
yy = ifelse(is.na(yy), 0, yy)
zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.1 & zero.p$zero.p<0.3, ])
yy = yy[, rownames(zero.p)]

##########################
## LMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
f1 = f2 = f3 = f4 = f5 = f6 = f7 = f8 = list()
out1 = out2 = out3 = out4 = out5 = out6 = out7 = out8 = out9 = out10 = out11 = out12 = matrix(NA, ncol(yy), 1)

for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f1 = lme(y0 ~ group+ weeks + age + race, random = ~ 1|subject)
    out1[j, ] = summary(f1)$tTable[2, 5]
    
    f2 = lme(y0 ~ group*weeks  + age + race, random = ~ 1|subject)
    out3[j, ] = summary(f2)$tTable[2, 5]
    out25[j, ] = summary(f2)$tTable[6, 5]
    
    f3 = lme(y0 ~ group*weeks + age + race, random = ~ 1|subject, correlation = corAR1())
    out5[j, ] = summary(f3)$tTable[2, 5]
    out7[j, ] = summary(f3)$tTable[6, 5]
    
    f4 = lme(y0 ~ group*weeks  + age + race, random = list(subject = pdDiag(~weeks)))
    out9[j, ] = summary(f4)$tTable[2, 5]
    out11[j, ] = summary(f4)$tTable[6, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
out_lmm <- cbind(out1, out3, out5, out7, out9, out11)

##########################
## ZIGMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
#di = matrix(NA, ncol(yy), 3)
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f5 = lme.zig(y1 ~ group + weeks  + age + race + offset(log(N))|1, family = "zig", random = ~ 1|subject)
    out2[j, ] = summary(f5)$tTable[2, 5]
    
    f6 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N))|1, family = "zig", random = ~ 1|subject)
    out4[j, ] = summary(f6)$tTable[2, 5]
    out26[j, ] = summary(f6)$tTable[6, 5]
    
    f7 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N))|1, family = "zig", random = ~ 1|subject, correlation = corAR1())
    out6[j, ] = summary(f7)$tTable[2, 5]
    out8[j, ] = summary(f7)$tTable[6, 5]
    
    f8 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N))|1, family = "zig", random = list(subject = pdDiag(~weeks)))
    out10[j, ] = summary(f8)$tTable[2, 5]
    out12[j, ] = summary(f8)$tTable[6, 5]

    f13 = lme.zig(y1 ~ group + weeks  + age + race + offset(log(N)), family = "zig", random = ~ 1|subject)
    out19[j, ] = summary(f13)$tTable[2, 5]
    
    f14 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N)), family = "zig", random = ~ 1|subject)
    out20[j, ] = summary(f14)$tTable[2, 5]
    out27[j, ] = summary(f14)$tTable[6, 5]
    
    f15 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N)), family = "zig", random = ~ 1|subject, correlation = corAR1())
    out21[j, ] = summary(f15)$tTable[2, 5]
    out22[j, ] = summary(f15)$tTable[6, 5]
    
    f16 = lme.zig(y1 ~ group*weeks  + age + race + offset(log(N)), family = "zig", random = list(subject = pdDiag(~weeks)))
    out23[j, ] = summary(f16)$tTable[2, 5]
    out24[j, ] = summary(f16)$tTable[6, 5]

    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
out_zigmm <- cbind(out2, out4, out6, out8, out10, out12)
```
