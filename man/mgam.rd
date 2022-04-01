
\name{mgam}
\Rdversion{1.0}
\alias{mgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting Mixed Models Separately for Many Responses Using \code{\link{gam}} in the package \bold{mgcv}
}

\description{
  This function fits mixed models separately for many responses using \code{\link{gam}} in the package \bold{mgcv}.    
}

\usage{  
mgam(y, formula, data, family=nb(), paraPen=NULL, min.p=0, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  a matrix of responses; each column is a response and rows are samples. 
  }
  \item{formula}{ 
  a one-sided formula of the form \code{~ x} (i.e., the respose is omitted); the right side of \code{~} is the same as in \code{\link{gam}}. 
}
  \item{data, family, paraPen}{ 
  the same as in \code{\link{gam}}.
}
\item{min.p}{
  a value in [0, 1). The responses with the proportion of non-zero values > min.p are analyzed.
}
\item{verbose}{
  logical. If \code{TRUE}, print out computational time.
}

}

\details{
 This function analyzes the responses in \code{y} by repeated calls to \code{\link{gam}} in the package \bold{mgcv}. 
 
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
  \code{\link{gam}}, \code{\link{mms}}, \code{\link{mglmmTMB}}  
}

\examples{
library(NBZIMM)

data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)

N = sam[, "Total.Read.Counts"]        
Days = sam$GA_Days
Age = sam$Age
Race = sam$Race
preg = sam$pregnant; table(preg)
subject = sam[, "Subect_ID"]; table(subject)

# analyze all taxa with a given nonzero proportion

f = mgam(y=otu, formula= ~ offset(log(N)) + Days + Age + Race + preg + s(subject, bs="re"), 
         family=nb(), min.p=0.2)
out = fixed(f)

para = out$parametric_terms
para = para[para[,2]!="(Intercept)", ]
par(mfrow = c(1, 1), cex.axis = 1, mar = c(2, 12, 4, 4))
plot.fixed(res=para[,c("Estimate","Std.Error","padj")], 
           threshold=0.001, gap=500, col.pts=c("black", "grey"),
           cex.axis=0.6, cex.var=0.7)

preg.coefs = para[para$variables=="preg",]
par(mfrow = c(1, 1), cex.axis = 1, mar = c(2, 10, 4, 4))
plot.fixed(res=preg.coefs[,c("Estimate","Std.Error","padj")], 
           threshold=0.05, gap=300, main="Covariate: pregnant",
           cex.axis=0.6, cex.var=0.7)

g = heat.p(df=para, p.breaks = c(0.001, 0.01, 0.05), 
       colors = c("black", "darkgrey", "grey", "lightgrey"),
       zigzag=c(T,F), abbrv=c(T,F), margin=c(1,6), y.size=8,
       legend=T)
}