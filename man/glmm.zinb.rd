
\name{glmm.zinb}
\Rdversion{1.0}
\alias{glmm.zinb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Zero-Inflated Negative Binomial Mixed Models
}

\description{
  This function sets up and fits zero-inflated negative binomial mixed models for analyzing zero-inflated count responses with multilevel data structures (for example, clustered data and longitudinal studies).    
}

\usage{  
glmm.zinb(fixed, random, data, correlation, 
        zi_fixed = ~1, zi_random = NULL,
        niter = 30, epsilon = 1e-05, verbose = TRUE, ...)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{ 
  a formula for the fixed-effects part of the negative binomial model, including the count outcome. This argument is the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
  \item{random}{ 
  a formula for the random-effects part of the negative binomial model. It only contain the right-hand side part, e.g., ~ time | id, where time is a variable, and id the grouping factor. This argument is the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
  \item{data}{ 
  a data.frame containing all the variables. 
}
  \item{correlation}{ 
  an optional correlation structure. It is the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
  \item{zi_fixed, zi_random}{
  formulas for the fixed and random effects of the zero inflated part. only contain the right-hand side part. 
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
  Zero-inflated models are two-component mixture models combining a point mass at zero with a negative binomial distribution for count response. Thus, there are two sources of zeros: zeros may come from either the point mass at zero or the negative binomial distribution. For modeling the unobserved state, a logistic model (or logistic mixed model) is used.
  
  The function is an alteration of the function \code{\link{glmmPQL}} in the package \bold{MASS}, which fits generalized linear mixed models using Penalized Quasi-Likelihood and works by repeated calls to the function \code{\link{lme}} in the package \bold{nlme}. It incorporates EM algorithms into the procedure of \code{\link{glmmPQL}} to fit zero-inflated mixed models for analyzing zero-inflated count responses in mutilevel study designs, for example, clustered and longitudinal studies. The function allows for multiple and correlated group-specific (random) effects (the argument \code{random}) and various types of within-group correlation structures (the argument \code{correlation}) described by \code{\link{corStruct}} in the package \bold{nlme}. 
}

\value{
  A fitted model object of class \code{lme} inheriting from \code{\link{lme}}, which can be summarized by functions in the package \bold{nlme} for the count distribution part.
  
  The object contains additional components for the zero-inflation part: 
\item{theta}{the estimate of the dispersion parameter;} 
\item{zero.prob}{the zero-state probabilities;}
\item{zero.indicator}{the conditional expectations of the zero indicators;}
\item{fit.zero}{the fitted logistic (mixed) model for the zero-inflation part;}
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
  \code{\link{lme}}, \code{\link{glmmPQL}}, \code{\link{glmm.nb}}, \code{\link{lme.zig}} 
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

y = otu[, 1]

data = data.frame(y=y, Days=Days, Age=Age, Race=Race, preg=preg, N=N, subject=subject)
f = glmm.zinb(fixed = y ~ Days + Age + Race + preg + offset(log(N)), 
              random = ~ 1 | subject, data = data) 
summary(f)
fixed(f)
summary(f$fit.zero)

library(GLMMadaptive) # compared with GLMMadaptive 
data = data.frame(y=y, Days=Days, Age=Age, Race=Race, preg=preg, N=N, subject=subject)
f0 = mixed_model(fixed = y ~ Days + Age + Race + preg + offset(log(N)),
                 random = ~ 1 | subject, zi_fixed = ~1,
                 data = data, family = zi.negative.binomial())
summary(f0)
library(glmmTMB) # compared with glmmTMB
f0 = glmmTMB(y ~ Days + Age + Race + preg + offset(log(N)) + (1 | subject), data = data, 
            family = nbinom2, zi = ~ 1)
summary(f0)

}