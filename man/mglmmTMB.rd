
\name{mglmmTMB}
\Rdversion{1.0}
\alias{mglmmTMB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting Mixed Models Separately for Many Responses Using \code{\link{glmmTMB}} in the package \bold{glmmTMB}
}

\description{
  This function fits mixed models separately for many responses using \code{\link{glmmTMB}} in the package \bold{glmmTMB}.    
}

\usage{  
mglmmTMB(y, formula, data, family = nbinom2(), 
        zi = ~0, disp = ~1,
        min.p = 0, verbose = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  a matrix of responses; each column is a response and rows are samples. 
  }
  \item{formula}{ 
  a one-sided formula of the form \code{~ x} (i.e., the respose is omitted); the right side of \code{~} is the same as in \code{\link{glmmTMB}}. 
}
  \item{data, family, zi, disp}{ 
  the same as in \code{\link{glmmTMB}}.
}
\item{min.p}{
  a value in [0, 1). The responses with the proportion of non-zero values > min.p are analyzed.
}
\item{verbose}{
  logical. If \code{TRUE}, print out computational time.
}

}

\details{
 This function analyzes the responses in \code{y} by repeated calls to \code{\link{glmmTMB}} in the package \bold{glmmTMB}. 
 
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
  \code{\link{glmmTMB}}, \code{\link{mms}}, \code{\link{mgam}}  
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
f1 = mglmmTMB(y = otu, formula = ~ Days + Age + Race + preg + offset(log(N))+ (1 | subject), 
        data = data, family = nbinom2, min.p = 0.2)

# zero-inflated negative binomial mixed model
f2 = mglmmTMB(y = otu, formula = ~ Days + Age + Race + preg + offset(log(N))+ (1 | subject), 
              data = data, family = nbinom2, zi = ~1,
              min.p = 0.2)

# display the results
f = f1
out = fixed(f)
out = out$cond
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