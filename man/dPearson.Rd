\name{dPearson}
\alias{Pearson}
\alias{dPearson}
\alias{pPearson}
\alias{qPearson}
\alias{rPearson}
\alias{sPearson}
\title{The Pearson product moment correlation coefficient.}
\usage{
dPearson(x, N, rho = 0, log = FALSE)

pPearson(q, N, rho = 0, lower.tail = TRUE, log.p = FALSE)

qPearson(p, N, rho = 0, lower.tail = TRUE, log.p = FALSE)

rPearson(n, N, rho = 0)

sPearson(N, rho = 0)
}
\arguments{
  \item{x,q}{vector of sample correlations}

  \item{p}{vector of probabilities}

  \item{rho}{vector of population correlations}

  \item{N}{vector of numbers of observations, \eqn{(N >
  3)}}

  \item{n}{number of values to generate. If n is a vector,
  length(n) values will be generated}

  \item{log,log.p}{logical vector; if TRUE, probabilities p
  are given as log(p)}

  \item{lower.tail}{logical vector; if TRUE (default),
  probabilities are \eqn{P[R <= r]}, otherwise, \eqn{P[R >
  r]}}
}
\value{
The output values conform to the output from other such
functions in \R. \code{dPearson()} gives the density,
\code{pPearson()} the distribution function and
\code{qPearson()} its inverse. \code{rPearson()} generates
random numbers. \code{sPearson()} produces a list
containing parameters corresponding to the arguments --
mean, median, mode, variance, sd, third cental moment,
fourth central moment, Pearson's skewness, skewness, and
kurtosis.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for the distribution of
Pearson's product moment correlation.
}
\examples{
pPearson(0.5, N=10)
pPearson(q=0.5, N=10, rho=0.3)
sPearson(N=10)
plot(function(x)dPearson(x,N=10,rho=0.7),-1,1)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\keyword{distribution}

