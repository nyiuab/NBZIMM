# NBZIMM: Negative Binomial and Zero-Inflated Mixed Models, with Applications to Microbiome Data Analysis

# Overview

This R package provides functions for setting up and fitting negative binomial mixed models and zero-inflated negative binomial and Gaussian models. These functions allow for mutiple and correlated group-specific (random) effects and various types of within-group correlation structures as described in the core package nlme, and return objects that can be summarized by functions in nlme. The methods can be used to analyze overdispersed and zero-inflated count or continuous responses with multilevel data structures (for example, clustered and longitudinal studies). 

Author: Nengjun Yi nyi@uab.edu; Maintainer: Nengjun Yi nyi@uab.edu

# Installation

library(remotes)

install_github("nyiuab/NBZIMM", force=T, build_vignettes=T)
