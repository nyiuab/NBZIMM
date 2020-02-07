
\name{nonzero}
\Rdversion{1.0}
\alias{nonzero}

\title{
Proportions of Non-Zero Values (Prevalence)
}

\description{
This function calculates the proportion of non-zero values (prevalence) for each response and plots the proportions.  
}

\usage{
nonzero(y, total, min.p=0, sort=TRUE, plot=FALSE)
}

\arguments{
  \item{y}{ 
  a matrix of responses; each column is a response and rows are samples. 
  }
  \item{total}{
  optional. It is total number (total reads in microbiome data). If not provided, the number is calculated as the column sum of \code{y}.
  }
  \item{min.p}{
  a value in [0, 1). The responses with the proportion of non-zero values > min.p are returned.
}
\item{sort}{
  sort by the nonzero proportions of the responses into decreasing order.
}
  \item{plot}{
  plot zero proportion and total number or not.
  }
}

\details{

}

\value{
Proportions of non-zero values (i.e. prevalence) for responses: \code{nonzero.p}, total numbers for all samples: \code{total}, mean and standard deviation of total numbers: \code{total_mean_sd}, filtered responses: \code{y.filter}, and plots of ordered proportions of zero values and total numbers.
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

data(Romero)
names(Romero)
otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)
N = sam[, "Total.Read.Counts"]  # total reads

non = nonzero(y=otu, total=N, plot=T)
non

}

