\name{dNormScore}
\alias{NormScore}
\alias{dNormScore}
\alias{normOrder}
\alias{pNormScore}
\alias{qNormScore}
\alias{rNormScore}
\alias{sNormScore}
\title{Normal Scores distribution.}
\usage{
dNormScore(x, c, N, U, log = FALSE)

pNormScore(q, c, N, U, lower.tail = TRUE, log.p = FALSE)

qNormScore(p, c, N, U, lower.tail = TRUE, log.p = FALSE)

rNormScore(n, c, N, U)

sNormScore(c, N, U)
}
\arguments{
  \item{x,q}{vector of non-negative quantities}

  \item{p}{vector of probabilities}

  \item{n}{number of values to generate. If n is a vector,
  length(n) values will be generated}

  \item{c}{vector number of treatments}

  \item{N}{vector total number of observations}

  \item{U}{vector sum of reciprocals of the number of the c
  sample sizes}

  \item{log,log.p}{logical vector; if TRUE, probabilities p
  are given as log(p)}

  \item{lower.tail}{logical vector; if TRUE (default),
  probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X >
  x]}}
}
\value{
The output values conform to the output from other such
functions in \R. \code{dNormScore()} gives the density,
\code{pNormScore()} the distribution function and
\code{qNormScore()} its inverse. \code{rNormScore()}
generates random numbers. \code{sNormScore()} produces a
list containing parameters corresponding to the arguments
-- mean, median, mode, variance, sd, third cental moment,
fourth central moment, Pearson's skewness, skewness, and
kurtosis. \code{normOrder()} gives the expected values of
the normal order statistics for a sample of size N.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for the normal scores test.
A function to calculate expected values of normal order
statistics is included.
}
\details{
This is the Kruskal-Wallis statistic with ranks replaced by
the expected values of normal order statistics. There are c
treatments with sample sizes \eqn{n_j, j=1 \dots c}{nj, j=1
\ldots c}. The total sample size is \eqn{N=\sum_1^c
n_j}{N=Sum nj}. The distribution depends on c, N, and U,
where \eqn{U=\sum_1^c (1/n_j)}{U=Sum (1/nj)}.

Let \eqn{e_N(k)}{eN(k)} be the expected value of the
\eqn{k_{th}}{kth} smallest observation in a sample of N
independent normal variates. Rank all observations
together, and let \eqn{R_{ij}}{Rij} denote the rank of
observation \eqn{X_{ij}}{Xij}, \eqn{i=1 \dots n_j}{i=1
\ldots nj} for treatment \eqn{j=1 \dots c}{j=1 \ldots c},
then the normal scores test statistic is

\deqn{x=(N-1)\frac{1}{\sum_{k=1}^{N} e_N(k)^2}
\sum_{j=1}^{c}\frac{S_j^2}{n_j}}{x=(N-1)[1/Sum(1 \ldots
N)(eN(k)^2)]Sum(1 \ldots c)[(Sj^2)/nj]}

where \eqn{S_j=\sum_{i=1}^{n_j}(e_N(R_{ij}))}{Sj=Sum(1
\ldots nj)(eN(Rij))}.

See Lu and Smith (1979) for a thorough discussion and some
exact tables for small r and n.  The calculations made here
use an incomplete beta approximation -- the same one used
for Kruskal-Wallis, only differing in the calculation of
the variance of the statistic.

The expected values of the normal order statistics use a
modification of M.Maechler's C version of the Fortran
algorithm given by Royston (1982). Spot checking the values
against Harter (1969) confirms the accuracy to 4 decimal
places as claimed by Royston.
}
\examples{
#Assuming three treatments, each with a sample size of 5

pNormScore(2, 3, 15, 0.6)
pNormScore(c(0.11, 1.5, 5.6), 3, 15, 0.6) ## approximately 5\% 50\% and 95\%
sNormScore(3, 15, 0.6)
plot(function(x)dNormScore(x, c=5, N=15, U=0.6), 0, 5)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\references{
Harter, H.L. (1969). \emph{Order statistics and their use
in testing and estimation, volume 2.} U.S. Supp. of Doc.

Lu, H.T. and Smith, P.J. (1979) Distribution of normal
scores statistic for nonparametric one-way analysis of
variance. \emph{Jour. Am Stat. Assoc.} \bold{74.} 715-722.

Royston, J.P. (1982). Expected normal order statistics
(exact and approximate) AS 177. \emph{Applied Statistics.}
\bold{31.} 161-165.
}
\keyword{distribution}

