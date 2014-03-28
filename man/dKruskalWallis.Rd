\name{dKruskalWallis}
\alias{Kruskal}
\alias{KruskalWallis}
\alias{dKruskalWallis}
\alias{pKruskalWallis}
\alias{qKruskalWallis}
\alias{rKruskalWallis}
\alias{sKruskalWallis}
\title{Kruskall-Wallis distribution.}
\usage{
dKruskalWallis(x, c, N, U, log = FALSE)

pKruskalWallis(q, c, N, U, lower.tail = TRUE, log.p = FALSE)

qKruskalWallis(p, c, N, U, lower.tail = TRUE, log.p = FALSE)

rKruskalWallis(n, c, N, U)

sKruskalWallis(c, N, U)
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
functions in \R. \code{dKruskalWallis()} gives the density,
\code{pKruskalWallis()} the distribution function and
\code{qKruskalWallis()} its inverse.
\code{rKruskalWallis()} generates random numbers.
\code{sKruskalWallis()} produces a list containing
parameters corresponding to the arguments -- mean, median,
mode, variance, sd, third cental moment, fourth central
moment, Pearson's skewness, skewness, and kurtosis.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for the Kruskal-Wallis test.
}
\details{
This is a one-way layout with, perhaps, unequal sample
sizes for each treatment. There are c treatments with
sample sizes \eqn{n_j, j=1 \dots c}{nj, j=1 \ldots c}. The
total sample size is \eqn{N=\sum_1^c n_j}{N=Sum(1 \ldots
c)nj}. The distribution depends on c, N, and U, where
\eqn{U=\sum_1^c (1/n_j)}{U=Sum(1 \ldots c)(1/nj)}.

Let \eqn{R_j}{Rj} be the sum of the ranks for treatment
\eqn{j (j=1\dots c)}{j (j=1 \ldots c)} then the
Kruskal-Wallis statistic is
\deqn{x=\frac{12}{N(N+1)}\sum_{j=1}^{c}\frac{R_j^2}{n_j} -
3(N+1)}{x=(12/(N (N-1))Sum(1 \ldots c)(Rj/nj) - 3 (N+1)}

This is asymptotically equivalent to a chi-squared variable
with c-1 degrees of freedom.

The original paper is Kruskal and Wallis (1952) with errata
appearing in Kruskal and Wallis (1953).  No attempt is made
to calculate exact values, rather an incomplete beta
approximation is used following Wallace (1959).
}
\examples{
# Assuming three treatments, each with a sample size of 5.
pKruskalWallis(1, 3, 15, 0.6)
pKruskalWallis(c(.1,1.5,5.7), 3, 15, 0.6) ## approximately 5\% 50\% and 95\%
sKruskalWallis(3, 15, 0.6)
plot(function(x)dKruskalWallis(x, 3, 15, 0.6),0.5,8)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\references{
Kruskal, W.H. and Wallis, W.A. (1952) Use of ranks in
one-criterion variance analysis. \emph{Jour. Am. Stat.
Assoc.} \bold{47.} 583-634

Kruskal, W.H. and Wallis, W.A. (1953) Errata to Use of
ranks in one-criterion variance analysis. \emph{Jour. Am.
Stat. Assoc.} \bold{48.} 907-911.

Wallace, D.L. (1959). Simplified beta-approximations to the
Kruskal-Wallis H test. \emph{Jour. Am. Stat. Assoc.}
\bold{54.} 225-230.
}
\keyword{distribution}

