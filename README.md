# NBZIMM: Negative Binomial and Zero-Inflated Mixed Models, with Applications to Microbiome Data Analysis

# Overview

This R package provides functions for setting up and fitting negative binomial mixed models and zero-inflated negative binomial and Gaussian models. These functions allow for mutiple and correlated group-specific (random) effects and various types of within-group correlation structures as described in the core package nlme, and return objects that can be summarized by functions in nlme. The methods can be used to analyze overdispersed and zero-inflated count or continuous responses with multilevel data structures (for example, clustered and longitudinal studies). 

Author: Nengjun Yi nyi@uab.edu; Maintainer: Nengjun Yi nyi@uab.edu

# Installation

 library(remotes)

 install_github("nyiuab/NBZIMM", force=T, build_vignettes=T)

# Methods

# MODELS AND ALGORITHMS
      
We describe our models and algorithms with a two-level design where samples are grouped in subjects. Assume that a microbiome/metagenomic study collects \textit{n} subjects and \textit{n\textsubscript{i}} samples for the \textit{i}-th subject. For each sample, we measure the counts for \textit{m} taxa (OTU, species, genus, \textit{etc.}), \textit{C\textsubscript{ijh}}; \textit{h} = 1, $\cdot$  $\cdot$ $\cdot$ , \textit{m}, the total sequence read \textit{T\textsubscript{ij}}, and some relevant covariates \textit{X\textsubscript{ij}}. The goal is to detect significant associations between taxa and covariates. \par

The function \textit{glmm.nb} in \pkg{NBZIMM} allows us to analyze the data for taxon \textit{h} with NBMMs:


$$C_{ijh} \sim NB(C_{ijh} |\mu_{ijh}, \theta_h),\ log(\mu_{ijh}) = log(T_{ij})+X_{ij}\beta_h+G_{ij}b_{ih},\ b_{ih} \sim N(0, \Psi_h)$$ 


where the dispersion \textit{$\theta$\textsubscript{h }}determines the over-dispersion, the offset log(\textit{T\textsubscript{ij}}) accounts for the varying total sequence reads, \textit{$\beta$\textsubscript{h}} is a vector of fixed effects, and \textit{b\textsubscript{ih}} are random effects. The inclusion of the random effects accounts for subject-specific effects and avoids biased inference on the fixed effects [@Pinheiro2000]. \textit{glmm.nb} iteratively approximates the NBMMs by a linear mixed model, which in turn is fitted using the function \textit{lme} in the package \pkg{nlme}. The dispersion \textit{$\theta$\textsubscript{h}} is then updated using Newton-Raphson algorithm as in the function \textit{glm.nb} of \pkg{MASS}. This framework allows us to incorporate the powerful features of \textit{lme} into NBMMs.\ \ \  \par

The function \textit{glmm.zinb} implements ZINBMMs that directly model true zeros and can be more efficient for analyzing taxa with excessive zeros than \textit{glmm.nb}: \par


$$C_{ijh} \sim  \left\{
\begin{array}
{rr}
0 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ &\ with\ probability\ p_{ij} \\
NB(C_{ijh}|\mu_{ijh}, \theta_h) & with\ probability\ 1-p_{ij} 
\end{array}
\right.$$


Here, the means \textit{$ \mu$\textsubscript{ijh}} are modeled as above, and the zero-inflation probabilities \textit{p\textsubscript{ijh}} are assumed to depend on some covariates via a logistic regression logit(\textit{p\textsubscript{ijh}}) = \textit{Z\textsubscript{ij}$ \alpha$\textsubscript{h}} or logistic mixed model logit(\textit{p\textsubscript{ijh}}) = \textit{Z\textsubscript{ij}$ \alpha$\textsubscript{h}} + \textit{G\textsubscript{ij}a\textsubscript{ih}}, where \textit{$ \alpha$\textsubscript{h}} is a vector of fixed effects and the random effects \textit{a\textsubscript{ih}} $\sim N(0, \Phi_{h})$. \textit{glmm.zinb} employs a fast and stable EM-IWLS algorithm to fit the ZINBMMs, also taking advantage of the standard function \textit{lme }for analyzing linear mixed models in the package \pkg{nlme}. \par

Rather than directly analyzing the observed counts, some methods analyze transformed count data [@Paulson2013], for example, \textit{y\textsubscript{ijh}} = log\textsubscript{2}(\textit{C\textsubscript{ijh}} + 1). Also, in some bioinformatics pipeline, such as MetaPhlAn, the whole metagenome shotgun sequencing data are processed and output in terms of relative abundance as proportion data. The transformation for proportion data is commonly chosen as arcsine square root transformation 
${{y}_{ijh}}=\text{arcsine}\left(\sqrt{{{C}_{ijh}}/{{T}_{ij}}} \right)$
. The function \textit{lme.zig} in \pkg{NBZIMM} implements ZIGMMs for the transformed count or proportion data:


$$y_{ij} \sim  \left\{
\begin{array}
{rr}
0 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ & with\ probability\ p_{ij} \\
N(y_{ijh}|\mu_{ijh}, \sigma^2) & with\ probability\ 1-p_{ij} 
\end{array}
\right.$$


where the means \textit{$ \mu$\textsubscript{ijh}} and the zero-inflation probabilities \textit{p\textsubscript{ijh}} are modeled as in \textit{glmm.zinb}. The function \textit{lme.zig} uses a fast and stable EM algorithm to fit ZIGMMs.
