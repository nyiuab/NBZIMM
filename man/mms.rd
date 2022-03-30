
\name{mms}
\Rdversion{1.0}
\alias{mms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting Mixed Models Separately for Many Responses Using \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, or \code{\link{lme.zig}} 
}

\description{
  This function fits mixed models separately for many responses using \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, or \code{\link{lme.zig}}.    
}

\usage{  
mms(y, fixed, random, data, method = c("nb", "lme", "zinb", "zig"),
    correlation, zi_fixed = ~1, zi_random = NULL,
    niter = 30, epsilon = 1e-05, min.p = 0, verbose = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  a matrix of responses; each column is a response and rows are samples. 
  }
  \item{fixed}{ 
  a one-sided formula of the form \code{~ x} (i.e., the respose is omitted); the right side of \code{~} is the same as in \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, or \code{\link{lme.zig}}. 
}
  \item{data, random, correlation, zi_fixed, zi_random, niter, epsilon}{ 
  These arguments are the same as in \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, or \code{\link{lme.zig}}.
}
\item{method}{
  character specification of the method:
  \code{"nb"}: glmm.nb; \code{"lme"}: lme;       
  \code{"zinb"}: glmm.zinb; \code{"zig"}: lme.zig. 
}
\item{min.p}{
  a value in [0, 1). The responses with the proportion of non-zero values > min.p are analyzed.
}
\item{verbose}{
  logical. If \code{TRUE}, print out computational time.
}

}

\details{
 This function analyzes the responses in \code{y} by repeated calls to \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, or \code{\link{lme.zig}}. 
 
}

\value{
  A list including \code{fit},  \code{responses}, \code{variables}, and \code{call}:

  \item{fit}{fitted models for all the responses;}
  \item{responses}{names of all the responses;}
  \item{variables}{names of all covariates;}
}

\references{
  
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{lme}}, \code{\link{glmm.nb}}, \code{\link{glmm.zinb}}, \code{\link{lme.zig}}  
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

# analyze all taxa with a given nonzero proportion

data = data.frame(Days=Days, Age=Age, Race=Race, preg=preg, N=N, subject=subject)

# negative binomial mixed model
f1 = mms(y = otu, fixed = ~ Days + Age + Race + preg + offset(log(N)), 
        random = ~ 1 | subject, data = data,
        min.p = 0.2, method = "nb")

# zero-inflated negative binomial mixed model        
f2 = mms(y = otu, fixed = ~ Days + Age + Race + preg + offset(log(N)), 
        random = ~ 1 | subject, zi_fixed = ~1, data = data,
        min.p = 0.2, method = "zinb")
        
# zero-inflated gaussian mixed model        
f3 = mms(y = log2(otu+1), fixed = ~ Days + Age + Race + preg + offset(log(N)), 
        random = ~ 1 | subject, zi_fixed = ~1, data = data,
        min.p = 0.2, method = "zig")

# display the results
f = f1
out = fixed(f)
out = out$dist
out = out[out[,2]!="(Intercept)", ]

par(mfrow = c(1, 1), cex.axis = 1, mar = c(2, 12, 4, 4))
plot.fixed(res=out[,c("Estimate","Std.Error","padj")], 
           threshold=0.001, gap=500, col.pts=c("black", "grey"),
           cex.axis=0.6, cex.var=0.7)

preg.coefs = out[out$variables=="preg",]
par(mfrow = c(1, 1), cex.axis = 1, mar = c(2, 10, 4, 4))
plot.fixed(res=preg.coefs[,c("Estimate","Std.Error","padj")], 
           threshold=0.05, gap=300, main="Covariate: pregnant",
           cex.axis=0.6, cex.var=0.7)

g = heat.p(df=out, p.breaks = c(0.001, 0.01, 0.05), 
       colors = c("black", "darkgrey", "grey", "lightgrey"),
       zigzag=c(T,F), abbrv=c(T,F), margin=c(2.5,0.5), y.size=8,
       legend=T)
}