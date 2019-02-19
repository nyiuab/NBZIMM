
\name{nonzero}
\Rdversion{1.0}
\alias{nonzero}

\title{
Proportions of Non-Zero Values
}

\description{
This simple function calculates the proportion of non-zero values for each response and plots the proportions.  
}

\usage{
nonzero(y, total, plot = FALSE)
}

\arguments{
  \item{y}{ 
  a matrix of responses; each column is a response and rows are samples. 
  }
  \item{total}{
  optional. It is total number (total reads in microbiome data). If not provided, the number is calculated as the column sum of \code{y}.
  }
  \item{plot}{
  plot zero proportion and total number or not.
  }
}

\details{

}

\value{
ordered proportions of non-zero values for responses: \code{nonzero.p}, total numbers for all samples: \code{total}, mean and standard deviation of total numbers: \code{total_mean_sd}, and plots of ordered proportions of zero values and total numbers.
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

non = nonzero(y = otu, total = N, plot = T)
non

}

