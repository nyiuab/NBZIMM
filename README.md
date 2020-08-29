# NBZIMM: Negative Binomial and Zero-Inflated Mixed Models, with Applications to Microbiome Data Analysis

# Overview

This R package provides functions for setting up and fitting negative binomial mixed models and zero-inflated negative binomial and Gaussian models. These functions allow for mutiple and correlated group-specific (random) effects and various types of within-group correlation structures as described in the core package nlme, and return objects that can be summarized by functions in nlme. The methods can be used to analyze overdispersed and zero-inflated count or continuous responses with multilevel data structures (for example, clustered and longitudinal studies). 

Author: Nengjun Yi nyi@uab.edu; Maintainer: Nengjun Yi nyi@uab.edu
**License**: GPL

# Installation

 library(remotes)

 install_github("nyiuab/NBZIMM", force=T, build_vignettes=T)

# Methods

There are three methods available to analyze microbiome data in NBZIMM. In all three methods, we separately analyze each microbiome taxon. 

 ## Negative Binomial mixed models (NBMMs)
 
![](https://github.com/nyiuab/NBZIMM/blob/master/images/nbmms.PNG?raw=true)

 ## Zero-inflated Negative Binomial mixed models (ZINBMMs)

![](https://github.com/nyiuab/NBZIMM/blob/master/images/zinbmms.PNG?raw=true)

 ## Zero-inflated Gaussian mixed models (ZIGMMs)
 
![](https://github.com/nyiuab/NBZIMM/blob/master/images/zigmms.PNG?raw=true)

# Tutorials

For a tutorial of using the function *glmm.nb* to analyze microbiome data with **NBMMs** please see: 
[NBMMs](https://github.com/nyiuab/NBZIMM/tree/master/tutorial/nbmms.md)

For a tutorial of using the function *glmm.zinb* to analyze microbiome data with **ZINBMMs** please see: 
[ZINBMMs](https://github.com/nyiuab/NBZIMM/tree/master/tutorial/zinbmms.md)

For a tutorial of using the function *lme.zig* to analyze microbiome data with **ZIGMMs** please see: 
[ZIGMMs](https://github.com/nyiuab/NBZIMM/tree/master/tutorial/zigmms.md)



