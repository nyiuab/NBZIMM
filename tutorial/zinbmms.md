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

library(BhGLM)
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

N = clinical[, "Total.Read.Countsb"]  # total reads
mean(N); sd(N)
mean(log(N)); sd(log(N))

subject = clinical[, "subject.ID"]
sample_id = clinical[, "sample.ID"]
CST = clinical[, "CST"]
clinical$group = group

## pregnant women vs non-pregnant women measurement time in weeks
days = clinical[, "calculated.GA..weeks."]

pheno = round(otu[, -c(1,2)]/100*N)

#pheno = pheno[rownames(clinical), ]
dim(pheno)
colnames(pheno)

yy = as.matrix(pheno)  
yy = ifelse(is.na(yy), 0, yy)
zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.2 & zero.p$zero.p<0.8, ])
yy = yy[, rownames(zero.p)]

##########################
## LMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
f1 = f2 = f3 = f4 = f5 = f6 = f7 = f8 = f9 = f10 = f11 = f12 = f13 = f14 = f15 = f16 = list()
f17 = f18 = f19 = f20 = f21 = f22 = f23 = f24 = list()
out1 = out2 = out3 = out4 = out5 = out6 = out7 = out8 = out9 = out10 = out11 = out12 = out13 = out14 = out15 = out16 = out17 = out18 = matrix(NA, ncol(yy), 1)
out19 = out20 = out21 = out22 = out23 = out24 = matrix(NA, ncol(yy), 1)
out25 = out26 = out27 = out28 = out29 = out30 = matrix(NA, ncol(yy), 1)
out31 = out32 = out33 = out34 = out35 = out36 = matrix(NA, ncol(yy), 1)
out37 = out38 = out39 = out40 = out41 = out42 = matrix(NA, ncol(yy), 1)

for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f1 = lme(y0 ~ group, random = ~ 1|subject)
    out1[j, ] = summary(f1)$tTable[2, 5]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f2 = lme(y0 ~ group*days, random = list(subject = pdDiag(~days)))
    out2[j, ] = summary(f2)$tTable[2, 5]
    out3[j, ] = summary(f2)$tTable[4, 5]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f3 = lme(y0 ~ group*days, random = ~ 1|subject)
    out4[j, ] = summary(f3)$tTable[2, 5]
    out5[j, ] = summary(f3)$tTable[4, 5]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f4 = lme(y0 ~ group*days, random = ~ 1|subject, correlation = corAR1()) #random = list(subject = pdDiag(~days)))
    out6[j, ] = summary(f4)$tTable[2, 5]
    out7[j, ] = summary(f4)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#out_lmm <- cbind(out1, out2, out4, out5, out6, out7)

for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    f9 = glmm.nb(y ~ group + offset(log(N)), random = ~ 1|subject)
    out8[j, ] = summary(f9)$tTable[2, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    f10 = glmm.nb(y ~ group*days + offset(log(N)), random = list(subject = pdDiag(~days)))
    out9[j, ] = summary(f10)$tTable[2, 5]
    out10[j, ] = summary(f10)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    f11 = glmm.nb(y ~ group*days + offset(log(N)), random = ~ 1|subject)
    out11[j, ] = summary(f11)$tTable[2, 5]
    out12[j, ] = summary(f11)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    f12 = glmm.nb(y ~ group*days+ offset(log(N)), random = ~ 1|subject, correlation = corAR1()) # random = list(subject = pdDiag(~days)))
    out13[j, ] = summary(f12)$tTable[2, 5]
    out14[j, ] = summary(f12)$tTable[4, 5]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#out_nbmm <- cbind(out13, out14, out15, out16, out17, out18)

##########################
## ZIGMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
#di = matrix(NA, ncol(yy), 3)
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f5 = lme.zig(y1 ~ group + offset(log(N))|1, random = ~ 1|subject)
    out15[j, ] = summary(f5)$tTable[2, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f6 = lme.zig(y1 ~ group*days + offset(log(N))|1, random = list(subject = pdDiag(~days)))
    out16[j, ] = summary(f6)$tTable[2, 5]
    out17[j, ] = summary(f6)$tTable[4, 5]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f7 = lme.zig(y1 ~ group*days + offset(log(N))|1, random = ~ 1|subject)
    out18[j, ] = summary(f7)$tTable[2, 5]
    out19[j, ] = summary(f7)$tTable[4, 5]    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f8 = lme.zig(y1 ~ group*days+ offset(log(N))|1, random = ~ 1|subject, correlation = corAR1()) #random = list(subject = pdDiag(~days)))
    out20[j, ] = summary(f8)$tTable[2, 5]
    out21[j, ] = summary(f8)$tTable[4, 5]    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f13 = lme.zig(y1 ~ group+ offset(log(N)), random = ~ 1|subject)
    out22[j, ] = summary(f13)$tTable[2, 5]    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f14 = lme.zig(y1 ~ group*days + offset(log(N)), random = list(subject = pdDiag(~days)))
    out23[j, ] = summary(f14)$tTable[2, 5]
    out24[j, ] = summary(f14)$tTable[4, 5]    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f15 = lme.zig(y1 ~ group*days + offset(log(N)), random = ~ 1|subject)
    out25[j, ] = summary(f15)$tTable[2, 5]
    out26[j, ] = summary(f15)$tTable[4, 5]
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    y1 = log2(y + 1)#asin(sqrt(y/N))
    f16 = lme.zig(y1 ~ group*days+ offset(log(N)), random = ~ 1|subject, correlation = corAR1()) #random = list(subject = pdDiag(~days)))
    out27[j, ] = summary(f16)$tTable[2, 5]
    out28[j, ] = summary(f16)$tTable[4, 5]    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#out_zigmm <- cbind(out2, out4, out6, out8, out10, out12)

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
