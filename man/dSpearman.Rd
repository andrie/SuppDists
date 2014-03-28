\name{dSpearman}
\alias{dSpearman}
\alias{pSpearman}
\alias{qSpearman}
\alias{rSpearman}
\alias{sSpearman}
\title{Spearman's rho.}
\usage{
dSpearman(x, r, log = FALSE)

pSpearman(q, r, lower.tail = TRUE, log.p = FALSE)

qSpearman(p, r, lower.tail = TRUE, log.p = FALSE)

rSpearman(n, r)

sSpearman(r)
}
\arguments{
  \item{x,q}{vector of non-negative quantities}

  \item{p}{vector of probabilities}

  \item{n}{number of values to generate. If n is a vector,
  length(n) values will be generated}

  \item{r}{(r >= 3) vector of number of observations}

  \item{log,log.p}{logical vector; if TRUE, probabilities p
  are given as log(p)}

  \item{lower.tail}{logical vector; if TRUE (default),
  probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X >
  x]}}
}
\value{
The output values conform to the output from other such
functions in \R. \code{dSpearman()} gives the density,
\code{pSpearman()} the distribution function and
\code{qSpearman()} its inverse. \code{rSpearman()}
generates random numbers. \code{sSpearman()} produces a
list containing parameters corresponding to the arguments
-- mean, median, mode, variance, sd, third cental moment,
fourth central moment, Pearson's skewness, skewness, and
kurtosis.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for Spearman's rho.
}
\details{
Spearman's rho is the rank correlation coefficient between
r pairs of items. It ranges from -1 to 1. Denote by d, the
sum of squares of the differences between the matched
ranks, then x is given by:

\deqn{1-\frac{6d}{r(r^2-1)}}{1-6 d /(r(r^2-1))}

This is, in fact, the product-moment correlation
coefficient of rank differences.  See Kendall (1975),
Chapter 2. It is identical to Friedman's chi-squared for
two treatments scaled to the -1, 1 range -- if X is the
Friedman statistic, then \eqn{\rho=frac{X}{r-1)-1}}{rho =
X/(r-1) -1}.

Exact calculations are made for \eqn{r \le 100}{r <= 100}

These exact calculations are made using the algorithm of
Kendall and Smith (1939).

The incomplete beta, with continuity correction, is used
for calculations outside this range.
}
\examples{
pSpearman(.95, 10)
pSpearman(c(-0.55,0,0.55), 10) ## approximately 5\% 50\% and 95\%
sSpearman(10)
plot(function(x)dSpearman(x, 10),-.9,.9)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\references{
Kendall, M. (1975). \emph{Rank Correlation Methods.}
Griffin, London.

Kendall, M. and Smith, B.B. (1939). The problem of m
rankings. \emph{Ann. Math. Stat.} \bold{10.} 275-287.
}
\keyword{distribution}

