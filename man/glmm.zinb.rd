
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
glmm.zinb(fixed, random, data, correlation, zi.random = FALSE, 
        niter = 30, epsilon = 1e-05, verbose = TRUE, ...)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{ 
  symbolic description of the fixed-effects part of the model, see details.
}
  \item{random, data, correlation}{ 
  These arguments are the same as in the function \code{\link{lme}} in the package \bold{nlme}.
}
\item{zi.random}{
  logical. If \code{TRUE}, include the random effect specified by \code{random} in the zero-inflation part. 
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
  Zero-inflated models are two-component mixture models combining a point mass at zero with a negative binomial distribution for count response. Thus, there are two sources of zeros: zeros may come from either the point mass at zero or the negative binomial distribution. For modeling the unobserved state, a logistic model (logistic mixed model if \code{zi.random = TRUE}) is used.
  
  The argument \code{fixed} is to specify the fixed-effects part in both components of the model: \code{fixed = y ~ x | z} giving the fixed-effects part of the distribution model y ~ x conditional on (|) the zero-inflation model ~ z. A different (or same) set of fixed-effects covariates could be used for the distribution component and zero-inflation component. If \code{fixed = y ~ x}, then the same covariates are employed in both components, equivalent to y ~ x | x. The simplest zero-inflation model only includes an intercept: y ~ x | 1, and thus all zeros have the same probability of belonging to the zero component. Offsets can be specified in both components of the model : y ~ x + offset(x1) | z + offset(z1).
  
  The function is an alteration of the function \code{\link{glmmPQL}} in the package \bold{MASS}, which fits generalized linear mixed models using Penalized Quasi-Likelihood and works by repeated calls to the function \code{\link{lme}} in the package \bold{nlme}. It incorporates EM algorithms into the procedure of \code{\link{glmmPQL}} to fit zero-inflated mixed models for analyzing zero-inflated count responses in mutilevel study designs, for example, clustered and longitudinal studies. The function allows for multiple and correlated group-specific (random) effects (the argument \code{random}) and various types of within-group correlation structures (the argument \code{correlation}) described by \code{\link{corStruct}} in the package \bold{nlme}. 
  
  Missing data may not be properly handled in some situations and thus should be removed before the analysis. 
}

\value{
  A fitted model object of class \code{lme} inheriting from \code{\link{lme}}, which can be summarized by functions in the package \bold{nlme} for the count distribution part.
  
  The object contains additional components for the zero-inflation part: 
\item{theta}{the estimate of the dispersion parameter;} 
\item{zero.prob}{the zero-state probabilities;}
\item{zero.indicator}{the conditional expectations of the zero indicators;}
\item{xz}{the design matrix of the zero-inflation part;}  
\item{offsetz}{the offset of the zero-inflation part;}
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

non = nonzero(y = otu, total = N, plot = F)
nonzero.p = non[[1]]

y = otu[, names(nonzero.p)[1]]

f = glmm.zinb(y ~ Days + Age + Race + preg + offset(log(N)) | 1, random = ~ 1 | subject) 
summary(f)
fixed(f)
summary(f$fit.zero)

}