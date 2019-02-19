
\name{heat.p}
\Rdversion{1.0}
\alias{heat.p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ggplot2 heatmap of p values  
}

\description{
Visualizing p values and positive coefficients for many responses and variables.
}

\usage{
heat.p(df, p.breaks = c(0.001, 0.01, 0.05), 
      colors = c("darkblue", "blue", "lightblue", "grey"),
      title = NULL, subtitle = NULL,
      filter = c(FALSE, FALSE), x.size, y.size, symbol.size,
      zigzag = c(FALSE, FALSE), abbrv = c(FALSE, FALSE), 
      margin = c(1, 1), legend = TRUE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
   data frame, giving response names, variable names, coefficient estimates and p values for many responses and variables. The column names must include "responses", "variables", "Estimate" and "pvalue". 
}
  \item{p.breaks}{
   a vector (not including 0 and 1), partitioning p values into intervals. 
}
\item{colors}{
  different colors for different p value intervals.
}
\item{title}{
  The text for title for the plot.
}
\item{subtitle}{
  The text for subtitle below the title.
}
\item{filter}{
 TRUE or FALSE. consist of two elements, the first is for responses (y-axis) and the second is for variables (x-axis). If \code{TRUE}, remove responses or variables with all p values larger than the maximun value of \code{p.breaks}.     
}
\item{x.size, y.size}{ 
  size for names of x-axis or y-axis.
}
\item{symbol.size}{ 
  size for the symbol "+" for positive estimates (only those with p values smaller than the maximun value of \code{p.breaks}).
}
\item{zigzag}{
 TRUE or FALSE. consist of two elements, the first is for responses (y-axis) and the second is for variables (x-axis). If \code{TRUE}, re-arrange y-axis or x-axis labels in a zigzag shape.     
}
\item{abbrv}{
 TRUE or FALSE. consist of two elements, the first is for responses (y-axis) and the second is for variables (x-axis). If \code{TRUE}, abbreviate y-axis or x-axis labels.     
}
\item{marigin}{
 a numeric vector of two values for the margins in centimeter. the first is for the left and the second is for the bottom.      
}
\item{legend}{
 TRUE or FALSE. show the legend for p values or not. 
}
  
}

\details{
  
}

\value{
 ggplot2 object
}

\author{
 Boyi Guo, boyiguo1@uab.edu; Xinyan Zhang, abbyyan3@gmail.com; Nengjun Yi, nyi@uab.edu 
}

\seealso{
  
}

\examples{
see examples in \code{\link{mms}} 
}

