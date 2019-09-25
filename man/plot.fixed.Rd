
\name{plot.fixed}
\Rdversion{1.0}
\alias{plot.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Estimates, Standard Deviations and P-values of Fixed Effects 
}

\description{
Plotting estimates, standard deviations and p-values of fixed effects.
}

\usage{
plot.fixed(res, threshold = 0.05, main = " ", col.pts = "black", gap = 0, 
          show.all.vars = FALSE, show.pvalues = TRUE, 
          cex.main = 0.9, xlim = NULL, cex.axis = 0.8, cex.var = 0.8, cex.pts = 1, 
          pch.pts = 20, type = "p", lwd = 1, lty = 1, line = 0,  
          add = FALSE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
   A matrix of the estimates, standard deviations and p-values of fixed effects for the response or variable. 
}
  \item{threshold}{
   a p-value to determine variables to be shown: variables with p-value < \code{threshold}. 
}
\item{col.pts}{
  color of points and segments, default is black. It can be a vector with two elements, which will use specified colors for different 
  variables determined by the augument \code{threshold}.
}
\item{gap}{
   a value for the distance between two flanking significant variables. 
}
\item{main, cex.main, xlim, type, lwd, lty, line}{ 
  These arguments are the same as in \code{\link[graphics]{plot}} and \code{\link[graphics]{par}}.
}
  \item{show.all.vars}{
  logical. If \code{TRUE},show names of all variables. The default is \code{FALSE}.
}
  \item{show.pvalues}{
  logical. If \code{TRUE}, show p-values for variables. The default is \code{TRUE}. The p-values are shown on the right side of the plot.
}
  \item{cex.var}{
  the fontsize of the varible names, default = \code{0.8}.
}
  \item{cex.pts}{
  the size of points, default = \code{1}.
}
  \item{pch.pts}{
  symbol of points, default is solid dot.
}
  \item{add}{
  logical. if \code{TRUE}, plot over the existing plot. The default is \code{FALSE}. 
}

}

\details{
This function plots the estimates of coefficients, intervals and p-values for a response or variable. It uses different colors to distinguish between significant and insignificant variables based on a threshold.  
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link[graphics]{plot}}, \code{\link[graphics]{par}}
}

\examples{
see examples in \code{\link{mms}} 
}

