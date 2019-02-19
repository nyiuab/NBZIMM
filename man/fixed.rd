
\name{fixed}
\Rdversion{1.0}
\alias{fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extracting and Summarizing Fixed Effects
}

\description{
This function is to extract the estimates, standard deviations and p-values of fixed effects  
from the output of \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, \code{\link{lme.zig}}, or \code{\link{mms}}.
}

\usage{
fixed(object) 

fixed.nb(object)

fixed.lme(object)

fixed.zi(object)

fixed.mms(object)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  an object from \code{\link{glmm.nb}}, \code{\link{lme}}, \code{\link{glmm.zinb}}, \code{\link{lme.zig}}, or \code{\link{mms}}. 
}
  
}

\details{

}

\value{
This function returns the following values for the distribution part (\code{dist}) and/or the zero-inflation part (\code{zero}): 

\item{Estimate}{intercept and coefficient estimates;} 
\item{Std. Error}{standard errors of estimates;}
\item{pvalue}{p values for testing coefficients;} 

}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{

}

\examples{
see examples in \code{\link{glmm.nb}}, \code{\link{glmm.zinb}}, \code{\link{lme.zig}}, and \code{\link{mms}} 
}

