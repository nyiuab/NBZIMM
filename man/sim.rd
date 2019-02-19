
\name{sim}
\Rdversion{1.0}
\alias{sim}

\title{
Simulating Overdispersed and/or Zero-Inflated Count Data with Random Effects
}

\description{
This function is to simulate overdispersed and/or zero-inflated count data with random effects.  
}

\usage{
sim(n.ind = 100, n.measure = 10, 
    x.d, coef.d = 1, tau.d = 1.5, 
    x.z, coef.z = 0, tau.z = 0, p.zero = 0, 
    theta = 2)
}

\arguments{
  \item{n.ind}{ 
  number of individuals (or clusters). 
  }
  \item{n.measure}{ 
  number of measures for each individual. 
  }
  \item{x.d}{ 
  design matrix of fixed-effects in the count part (distribution part). If not given, \code{x.d} will be generated as a binary predictor internally.
  }
  \item{coef.d}{ 
  coefficients of \code{x.d}.   
  }
  \item{tau.d}{ 
  standard deviation of random effect in the count part.   
  }
  \item{x.z}{ 
  design matrix of fixed-effects in the zero-inflation part.
  }
  \item{coef.z}{ 
  coefficients of \code{x.z}.   
  }
  \item{tau.z}{ 
  standard deviation of random effect in the zero-inflation part.   
  }
  \item{p.zero}{ 
  proportion of zeros from the the point mass at zero (i.e., not from the negative binomial count distribution). 
  }
  \item{theta}{ 
  shape parameter for negative binomial model of the count part.   
  }
  
}

\details{

}

\value{
a list containing individual ID \code{ind.ID}, design matrices \code{x.d} and \code{x.z}, Total number \code{T}, and count response \code{y}. 
}
\references{

}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  
}

\examples{

library(NBZIMM)

d = sim(n.ind = 100, n.measure = 10, coef.d = 1, tau.d = 1.5, 
        theta = 2, p.zero = 0)
ind = d$ind.ID
x = d$x.d
y = d$y
off = log(d$T)

f1 = glmm.nb(y ~ offset(off) + x, random = ~ 1|ind) 
summary(f1)
f1$theta

d = sim(n.ind = 100, n.measure = 10, coef.d = 1, tau.d = 1.5, 
        theta = 2, p.zero = 0.4)
ind = d$ind.ID
x = d$x.d
y = d$y
off = log(d$T)

f1 = glmm.zinb(y ~ offset(off) + x | 1, random = ~ 1|ind)
summary(f1)
f1$theta
unique(f1$zero.prob)

}

