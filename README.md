# NBZIMM: Negative Binomial and Zero-Inflated Mixed Models, with Applications to Microbiome Data Analysis

# Overview

This R package provides functions for setting up and fitting negative binomial mixed models and zero-inflated negative binomial and Gaussian models. These functions allow for mutiple and correlated group-specific (random) effects and various types of within-group correlation structures as described in the core package nlme, and return objects that can be summarized by functions in nlme. The methods can be used to analyze overdispersed and zero-inflated count or continuous responses with multilevel data structures (for example, clustered and longitudinal studies). 

Author: Nengjun Yi nyi@uab.edu; Maintainer: Nengjun Yi nyi@uab.edu

# Installation

Three ways to install the package in R:

1. Without Vignettes
```{r}
if (!requireNamespace("devtools")) install.packages("devtools")
library(devtools)
install_github("nyiuab/NBZIMM")
```
2. With Vignettes
```{r}
if (!requireNamespace("devtools")) install.packages("devtools")
if(!requireNamespace(R.rsp)) install.packages("R.rsp")
library(devtools)
install_github("nyiuab/NBZIMM", build_opts = c("--no-resave-data", "--no-manual"), force = T)
```
3. Download the NBZIMM zip file to your computer, and then install it to R.
