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

## Simulation Studies
The following R code is used for simulation studies (Setting 1) in a manuscript and the citation will added later.

```r
rm(list=ls(all=TRUE))
ls()

library(NBZIMM)
library(BhGLM)
library(nlme)
library(MASS)

### set number of individuals
n_1 <- 50
n_2 <- 100
n_3 <- 150 #100 

n <- n_1

### average intercept and slope
b_2 = c(0.2, 0.35)
b_3 = c(0.3, 0.8) 
b_1 = c(0, 0)

b0 = b_1
b1 = b_2

##### number of obserations parameters
m = 5
n.t = n*m
##### covariate matrix
X = NA

f1 = f2 = f3 = f4 = list()

bb = NULL

n.sims = 100
out1 = out2 = out3 = out4 = out5 = matrix(NA, n.sims, 3)
out6 = out7 = out8 = out9 = out10 = matrix(NA, n.sims, 3)

start.time <- Sys.time()

for (i in 1:n.sims){
  tryCatch({
    #set.seed(12345+i)
    diet = c(rep(1, n/2*m), rep(0, n/2*m))   
    id = c(rep(1:(n/2),m),rep((n/2+1):n,m))
    U_1 = rep(NA, n)
    sigma_1 = runif(1, 0.5, 1); #v0_1 = c(v0_1, sigma_1^2); sigma_1^2
    for (k in 1:n) U_1[k] = rnorm(1, 0, sigma_1)
    
    t = rep(log(1:m), each = n/2, times = 2)
    
    beta_1 = runif(1, b0[1], b0[2])
    beta_2 = runif(1, b0[1], b0[2])
    beta_3 = runif(1, b1[1], b1[2])
    
    inter = t*diet
    
    b.rep <- U_1[id]

    logT = runif(n.t, 7.1, 10.5)
    N = exp(logT)
    u0 = -7
    expmu0 = logT + u0 
    mu0 = log(expmu0)
    
    log.u <- mu0 + diet*beta_1 + t*beta_2 + beta_3*inter + b.rep# #+b.rep#
    u <- exp(log.u)
    
    nb.theta = runif(1, 0.1, 5)
    Y1 = rnegbin(n.t, mu = u, theta = nb.theta)
    Y = Y1
    
    X <- as.matrix(data.frame(id=id, t=t, diet=diet, y0=Y, N=N))
    
    ### set up data frame
    dat <- data.frame(X)
    
    N = dat$N
    y0 = dat$y0
    diet = dat$diet
    subject.ind = dat$id
    time.ind = dat$t
    
    y1 = asin(sqrt((y0)/(N)))
    
    f1 = lme(y1 ~ diet*time.ind,  random = ~1|subject.ind) 
    out1[i, ] = summary(f1)$tTable["diet:time.ind", ][c(1,2,5)]
    out6[i, ] = summary(f1)$tTable["diet", ][c(1,2,5)]
    
    f4 = glmm.nb(y0 ~ diet*time.ind + offset(log(N)), random = ~1|subject.ind, verbose = F) 
    out5[i, ] = summary(f4)$tTable["diet:time.ind", ][c(1,2,5)]
    out10[i, ] = summary(f4)$tTable["diet", ][c(1,2,5)]
       
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}
stop.time <- Sys.time()
round(difftime(stop.time, start.time, units = "min"), 3)

```

## Real Data Analysis
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
clinical <- clinical[clinical$Outcome != "Marginal", ] ### Same as in the paper
clinical <- clinical[clinical$V1 == "4", ] ### Same as in the paper

clinical <- clinical[clinical$Preg == "TRUE", ]
clinical <- clinical[clinical$BodySite == "Vaginal_Swab", ] ### To compare with CST analysis in the paper
otu <- otu[otu$Sample_ID %in% clinical$SampleID, ]

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
SubjectID = as.factor(clinical$SubjectID)
outcome = clinical$Outcome # preterm vs term groups
table(outcome)
term_label = ifelse(outcome == "Term", 0, 1)

weeks = as.numeric(log(clinical$GDColl)) # Days for samples collected
summary(weeks)

N = clinical$NumReads  # total reads
mean(N); sd(N)
mean(log(N)); sd(log(N))

##########################
## LMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
f1 = f2 = f3 = f4 = f5 = f6 = f7 = f8 = list()
out1 = out2 = out3 = out4 = out5 = out6 = out7 = out8 = out9 = out10 = out11 = out12 = matrix(NA, ncol(yy), 1)

for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f1 = lme(y0 ~ term_label, random = ~ 1|SubjectID)
    out1[j, ] = summary(f1)$tTable[2, 5]
    
    f2 = lme(y0 ~ term_label, random = list(SubjectID = pdDiag(~weeks)))
    out3[j, ] = summary(f2)$tTable[2, 5]
    
    f3 = lme(y0 ~ term_label*weeks, random = ~ 1|SubjectID)
    out5[j, ] = summary(f3)$tTable[2, 5]
    out7[j, ] = summary(f3)$tTable[4, 5]
    
    f4 = lme(y0 ~ term_label*weeks, random = list(SubjectID = pdDiag(~weeks)))
    out9[j, ] = summary(f4)$tTable[2, 5]
    out11[j, ] = summary(f4)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
out_lmm <- cbind(out1, out3, out5, out7, out9, out11)

##########################
## NBMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    #y0 = log2(y + 1)#asin(sqrt(y/N))
    f5 = glmm.nb(y ~ term_label + offset(log(N)), random = ~ 1|SubjectID)
    out2[j, ] = summary(f5)$tTable[2, 5]
    
    f6 = glmm.nb(y ~ term_label + offset(log(N)), random = list(SubjectID = pdDiag(~weeks)))
    out4[j, ] = summary(f6)$tTable[2, 5]
    
    f7 = glmm.nb(y ~ term_label*weeks + offset(log(N)), random = ~ 1|SubjectID)
    out6[j, ] = summary(f7)$tTable[2, 5]
    out8[j, ] = summary(f7)$tTable[4, 5]
    
    f8 = glmm.nb(y ~ term_label*weeks + offset(log(N)), random = list(SubjectID = pdDiag(~weeks)))
    out10[j, ] = summary(f8)$tTable[2, 5]
    out12[j, ] = summary(f8)$tTable[4, 5]
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
out_nbmm <- cbind(out2, out4, out6, out8, out10, out12)

out <- cbind(out1, out2, out3,  out4, out5, out6, out7, out8, out9, out10, out11, out12)
rownames(out) = colnames(yy)
est = out2[complete.cases(out2), 1]
est = est
res = sim.out(coefs.p = t(out), coefs.est = t(est), alpha = c(0.05, 0.01, 0.005, 0.001))
res
res2 = sim.out(coefs.p = t(p.adjust), coefs.est = t(est), alpha = c(0.05, 0.01, 0.005, 0.001))
res2
```

