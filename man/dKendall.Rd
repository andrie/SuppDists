\name{dKendall}
\alias{dKendall}
\alias{pKendall}
\alias{qKendall}
\alias{rKendall}
\alias{sKendall}
\title{The distribution of Kendall's tau.}
\usage{
dKendall(x, N, log = FALSE)

pKendall(q, N, lower.tail = TRUE, log.p = FALSE)

qKendall(p, N, lower.tail = TRUE, log.p = FALSE)

rKendall(n, N)

sKendall(N)
}
\arguments{
  \item{x,q}{vector of non-negative quantities}

  \item{p}{vector of probabilities}

  \item{n}{number of values to generate. If n is a vector,
  length(n) values will be generated}

  \item{N}{vector number of treatments}

  \item{log,log.p}{logical vector; if TRUE, probabilities p
  are given as log(p)}

  \item{lower.tail}{logical vector; if TRUE (default),
  probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X >
  x]}}
}
\value{
The output values conform to the output from other such
functions in \code{dKendall()} gives the density,
\code{pKendall()} the distribution function and
\code{qKendall()} its inverse. \code{rKendall()} generates
random numbers. \code{sKendall()} produces a list
containing parameters corresponding to the arguments --
mean, median, mode, variance, sd, third cental moment,
fourth central moment, Pearson's skewness, skewness, and
kurtosis.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for Kendall's tau.
}
\details{
There are two categories with N treatments each.  The
treatments are ranked for each category, and then sorted
according to the ranks for the first category.  This
produces a 2 by N array in which the numbers in the first
row are increasing from 1 to N.  The array is scanned, and
every time two adjacent ranks in the second row are not in
order, they are exchanged.  The scanning is repeated until
the second row is in increasing order.  Let s denote the
number of exchanges, then Kendall's tau is given by

\deqn{\tau=1-\frac{4s}{N(N-1)}}{tau=1-4 s/(N (N-1))}

This too is a product-moment correlation coefficient.  See
Kendall (1975), Chapter 2.  Other methods for calculating
the statistic are also discussed there.

The calculated values are exact for \eqn{N < 13},
thereafter an Edgeworth expansion is used.
}
\examples{
pKendall(0, N=10)
pKendall(c(-.42, 0.02, .42), N=10) ## approximately 5\% 50\% and 95\%
qKendall(.95, N=c(10, 20))
sKendall(N=10)
plot(function(x) dKendall(x, N=10), -0.5, 0.5)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\references{
Kendall, M. (1975). \emph{Rank Correlation Methods.}
Griffin, London.
}
\keyword{distribution}

