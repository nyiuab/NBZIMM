
\name{get.fixed}
\Rdversion{1.0}
\alias{get.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract Fixed Effects for A Response or Variable from a 'mms' Output  
}

\description{
  This function extracts the estimates, standard deviations and p-values of fixed effects for a response or variable from the output of \code{\link{mms}}.    
}

\usage{  
get.fixed(object, part=c("dist", "zero"), vr.name, sort=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object from \code{\link{mms}}. 
}
  \item{part}{
  distribution part or zero-inflation part. 
}
  \item{vr.name}{
  name of a variable or response. 
}
  \item{sort}{
  sort by the adjusted p-values into ascending order.
  }

}

\details{
 
}

\value{
  A matrix consists of the estimates, standard deviations, p-values and adjusted p-values of fixed effects for the response or variable.
}

\references{
  
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{mms}} 
}

\examples{
see examples in \code{\link{mms}} 

}