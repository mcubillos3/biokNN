% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extra.R
\name{create_multilevel}
\alias{create_multilevel}
\title{Generate multilevel dataset}
\usage{
create_multilevel(
  nClass = 10,
  nVars = 1,
  classMean = 10,
  classSD = 0,
  beta0 = 0,
  tau0 = 1,
  beta = c(1),
  tau = c(1),
  sigma2 = 1
)
}
\arguments{
\item{nClass}{number of classes}

\item{nVars}{number of independent variables (X)}

\item{classMean}{average number of observations per class}

\item{classSD}{standard deviation of the number of observations per class}

\item{beta0}{intercept parameter}

\item{tau0}{variance of the parameter between classes}

\item{beta}{vector with the slope parameters, one for each independent variable}

\item{tau}{vector with the variance of the slope parameters, one for each independent variable}

\item{sigma2}{error variance}
}
\value{
A dataframe with the multilevel dataset
}
\description{
This function returns a dataframe with a multilevel structure. It generates a dataframe using a varying
intercepts/varying slopes linear regression with a single target variable y.
}
\examples{
df <- create_multilevel(nClass = 20,
                          nVars = 1,
                          classMean = 10,
                          classSD = 2)
}
