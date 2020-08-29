# NBZIMM: Negative Binomial and Zero-Inflated Mixed Models, with Applications to Microbiome Data Analysis

# Overview

This R package provides functions for setting up and fitting negative binomial mixed models and zero-inflated negative binomial and Gaussian models. These functions allow for mutiple and correlated group-specific (random) effects and various types of within-group correlation structures as described in the core package nlme, and return objects that can be summarized by functions in nlme. The methods can be used to analyze overdispersed and zero-inflated count or continuous responses with multilevel data structures (for example, clustered and longitudinal studies). 

Author: Nengjun Yi nyi@uab.edu; Maintainer: Nengjun Yi nyi@uab.edu

# Installation

 library(remotes)

 install_github("nyiuab/NBZIMM", force=T, build_vignettes=T)

# Methods

There are three methods available to analyze microbiome data in NBZIMM. In all three methods, we separately analyze each microbiome taxon. 
1. Negative Binomial mixed models (NBMMs)
<img src="https://github.com/nyiuab/NBZIMM/tree/master/pics/nbmms.PNG" width="600" align="center">
![GitHub Logo](/images/nbmms.png)
2. ZINBMMs
3. ZIGMMs


# Tutorials

[NBMMs](https://github.com/nyiuab/NBZIMM/tree/master/tutorial/nbmms.md)

[ZINBMMs](https://github.com/nyiuab/NBZIMM/tree/master/tutorial/zinbmms.md)

[ZIGMMs](https://github.com/nyiuab/NBZIMM/tree/master/tutorial/zigmms.md)



**License**: GPL
