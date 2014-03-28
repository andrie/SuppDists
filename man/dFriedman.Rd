\name{dFriedman}
\alias{dFriedman}
\alias{pFriedman}
\alias{qFriedman}
\alias{rFriedman}
\alias{sFriedman}
\title{Friedman's chi-square}
\usage{
dFriedman(x, r, N, log = FALSE)

pFriedman(q, r, N, lower.tail = TRUE, log.p = FALSE)

qFriedman(p, r, N, lower.tail = TRUE, log.p = FALSE)

rFriedman(n, r, N)

sFriedman(r, N)
}
\arguments{
  \item{x,q}{vector of non-negative quantities}

  \item{p}{vector of probabilities}

  \item{n}{number of values to generate. If n is a vector,
  \code{length(n)} values will be generated}

  \item{r}{vector of number of treatments}

  \item{N}{(N >= 2) vector of number of replications of
  each treatment}

  \item{log,log.p}{logical vector; if TRUE, probabilities p
  are given as log(p)}

  \item{lower.tail}{logical vector; if TRUE (default),
  probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X >
  x]}}
}
\value{
The output values conform to the output from other such
functions in \code{dFriedman()} gives the density,
\code{pFriedman()} the distribution function and
\code{qFriedman()} its inverse. \code{rFriedman()}
generates random numbers. \code{sFriedman()} produces a
list containing parameters corresponding to the arguments
-- mean, median, mode, variance, sd, third cental moment,
fourth central moment, Pearson's skewness, skewness, and
kurtosis.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for Friedman's chi square.
}
\details{
The Freidman chi-squared is used for nonparametric ANOVA.
The data in N rows of an \eqn{N \times r}{N x r} table are
ranked separately such that the ranks take the values from
1 to r in the N different rows. The distributions are
obtained on the assumption that there is no relationship
between the N rows.

Formulae:

Let \eqn{R_j}{Rj} be the sum of ranks for treatment \eqn{j
(j=1\dots r)}{j (j=1 \ldots r)}, then the Friedman
statistic is

\deqn{ }{x=[12/(N r (r+1))]Sum(1 \ldots r)(Rj^2) - 3 N
(r+1)}\deqn{x=\frac{12}{N r (r+1)}\sum_{j=1}^{r}R_j^2
-3N(r+1)}{x=[12/(N r (r+1))]Sum(1 \ldots r)(Rj^2) - 3 N
(r+1)}

this is asymptotically equivalent to a \eqn{\chi^2}{chi
squared} random variable. One may also calculate the chi
squared statistic for the usual analysis of variance which
gives

\deqn{ }{F=((N-1) x)/(N
(r-1)-x)}\deqn{F=\frac{(N-1)x}{N(r-1)-x}}{F=((N-1) x)/(N
(r-1)-x)}

which may be used with the F distribution functions in R
for degrees of freedom \eqn{(r-1)} and \eqn{(N-1)(r-1)}.
}
\note{
Exact calculations are made for the following values of the
parameters:

\tabular{ll}{ r \tab N \cr 2 \tab 100 \cr 3 \tab 30 \cr 4
\tab 15 \cr 5 \tab 8 } These exact calculations are made
using the algorithm of Kendall and Smith (1939).

The incomplete beta, with continuity correction, is used
for calculations outside these ranges.  Some appreciation
for the accuracy of the approximation may be obtained by
comparing the calculated values with exact tables such as
Odeh (1977).  Iman and Davenport (1980) have studied the
accuracy of various approximations.
}
\examples{
pFriedman(2, r=5, N=10)
pFriedman(c(.8,3.5,9.3), r=5, N=10) ## approximately 5\% 50\% and 95\%
sFriedman(r=5, N=10)
plot(function(x)dFriedman(x, r=5, N=10),0,10)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\references{
Kendall, M. and Smith, B.B. (1939). The problem of m
rankings. \emph{Ann. Math. Stat.} \bold{10.} 275-287.

Iman, R.L. and Davenport, J.M. (1980). Approximations of
the critical region of the Friedman statistic. \emph{Comm.
Stat. Theor. Meth.} \bold{A9(6).} 571-595.

Odeh, R.E. (1977). Extended tables of the distribution of
Friedman's S-statistic in the two-way layout. \emph{Commun.
Statist.-Simula. Computa.} \bold{B6(1).} 29-48.
}
\keyword{distribution}

