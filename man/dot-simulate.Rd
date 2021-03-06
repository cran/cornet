% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{.simulate}
\alias{.simulate}
\title{Data simulation}
\usage{
.simulate(n, p, cor = 0, prob = 0.1, sd = 1, exp = 1, frac = 1)
}
\arguments{
\item{n}{sample size\strong{:}
positive integer}

\item{p}{covariate space\strong{:}
positive integer}

\item{cor}{correlation coefficient \strong{:}
numeric between \eqn{0} and \eqn{1}}

\item{prob}{effect proportion\strong{:}
numeric between \eqn{0} and \eqn{1}}

\item{sd}{standard deviation\strong{:}
positive numeric}

\item{exp}{exponent\strong{:}
positive numeric}

\item{frac}{class proportion\strong{:}
numeric between \eqn{0} and \eqn{1}}
}
\value{
Returns invisible list with elements \code{y} and \code{X}.
}
\description{
Simulates data for unit tests
}
\details{
For simulating correlated features (\code{cor}\eqn{>0}),
this function requires the R package MASS
(see \code{\link[MASS]{mvrnorm}}).
}
\examples{
data <- cornet:::.simulate(n=10,p=20)
names(data)

}
