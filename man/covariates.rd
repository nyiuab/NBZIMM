
\name{covariates}
\Rdversion{1.0}
\alias{covariates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforming Covariates and Filling in Missing Values
}

\description{
  This function standardizes continuous variables, transforms categorical variables to indicator variables and centers the transformed variables.     
}

\usage{
covariates(x.con, x.cat, con.rescale = TRUE, cat.center = FALSE, fill.missing = TRUE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.con}{ 
  a matrix or data frame consisting of continuous covariates.    
}
  \item{x.cat}{ 
  a matrix or data frame consisting of categorical covariates.    
}
  \item{con.rescale}{ 
  if \code{con.rescale = TRUE}, standardize continuous variables. 
}
  \item{cat.center}{ 
  if \code{cat.center = TRUE}, center categorical variables. 
}
  \item{fill.missing}{
  If \code{fill.missing = TRUE}, missing values of a variable are replaced by the mean of the observed data.
  }
}

\details{

}

\note{

}

\value{
This function returns a data frame containing the rescaled variables.         
}

\references{

}

\author{
Nengjun Yi, nyi@uab.edu
}

\examples{
# fake data
age = rnorm(50, 30, 0.1)
sex = sample(x = c("male", "female"), size = 50, replace = T, prob = c(0.5, 0.5))
diet = sample(x = c(1, 5), size = 50, replace = T, prob = c(0.3, 0.7))
race = sample(c("Asian", "White", "Black"), size = 50, replace = T, prob = c(0.2, 0.4, 0.4))

x.con = age # continuous variable
x.cat = cbind(sex, diet, race) # categorical variables

x.con[1:2] = NA
x.cat[1,] = NA

x1 = covariates(x.con = x.con, x.cat = x.cat, con.rescale = F, cat.center = F, 
                fill.missing = T)
x1
x2 = covariates(x.con = x.con, x.cat = x.cat, con.rescale = T, cat.center = F, 
                fill.missing = T)
x2
x3 = covariates(x.con = x.con, x.cat = x.cat, con.rescale = T, cat.center = T, 
                fill.missing = T)
x3
}
