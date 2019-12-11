
\name{glmm.nb}
\Rdversion{1.0}
\alias{glmm.nb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Negative Binomial Mixed Models
}

\description{
  This function sets up and fits negative binomial mixed models for analyzing overdispersed count responses with multilevel data structures (for example, clustered data and longitudinal studies).    
}

\usage{  
glmm.nb(fixed, random, data, correlation,  
        niter = 30, epsilon = 1e-05, verbose = TRUE, ...)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{ 
  a formula for the fixed-effects part, including the outcome. This argument is the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
  \item{random}{ 
  a formula for the random-effects part. It only contain the right-hand side part, e.g., ~ time | id, where time is a variable, and id the grouping factor. This argument is the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
  \item{data}{ 
  a data.frame containing all the variables. 
}
  \item{correlation}{ 
  an optional correlation structure. It is the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
  \item{niter}{
  maximum number of iterations. 
}
 \item{epsilon}{
  positive convergence tolerance.
}
  \item{verbose}{
  logical. If \code{TRUE}, print out number of iterations and computational time.
}
  \item{...}{
  further arguments for \code{\link{lme}}.  
}

}

\details{
  This function is an alteration of the function \code{\link{glmmPQL}} in the package \bold{MASS}, which fits generalized linear mixed models using Penalized Quasi-Likelihood and works by repeated calls to the function \code{\link{lme}} in the package \bold{nlme}. It sets up and fits negative binomial mixed-effects models, which are the standard models for analyzing overdispersed count responses in mutilevel study designs, for example, clustered and longitudinal studies. The function allows for multiple and correlated group-specific (random) effects (the argument \code{random}) and various types of within-group correlation structures (the argument \code{correlation}) described by \code{\link{corStruct}} in the package \bold{nlme}.    
}

\value{
  A fitted model object of class \code{lme} inheriting from \code{\link{lme}}, which can be summarized by functions in the package \bold{nlme}. The object also contains the estimate of the dispersion parameter \code{theta}.
}
\references{
  Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer. 

  Venables, W. N. and Ripley, B. D. (2002) "Modern Applied Statistics with S". Fourth edition. Springer. 

  Xinyan Zhang, Himel Mallick, Xiangqin Cui, Andrew K. Benson, and Nengjun Yi (2017) Negative Binomial Mixed Models for Analyzing Microbiome Count Data. BMC Bioinformatics 18(1):4.
  
  Xinyan Zhang, Yu-Fang Pei, Lei Zhang, Boyi Guo, Amanda Pendegraft, Wenzhuo Zhuang and Nengjun Yi (2018) Negative Binomial Mixed Models for Analyzing Longitudinal Microbiome Data. Frontiers in Microbiology.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{lme}}, \code{\link{glmmPQL}} 
}

\examples{

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

non = nonzero(y = otu, total = N, plot = F)
nonzero.p = non[[1]]

y = otu[, names(nonzero.p)[1]]

f = glmm.nb(y ~ Days + Age + Race + preg + offset(log(N)), random = ~ 1 | subject)
summary(f)
fixed(f)


library(lme4)  # compared with lme4
f0 = glmer.nb(y ~ Days + Age + Race + preg + offset(log(N)) + (1 | subject))   
summary(f0)

}