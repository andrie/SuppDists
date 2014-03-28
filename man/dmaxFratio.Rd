\name{dmaxFratio}
\alias{dmaxFratio}
\alias{maxFratio}
\alias{pmaxFratio}
\alias{qmaxFratio}
\alias{rmaxFratio}
\alias{smaxFratio}
\title{The maximum F-ratio distribution.}
\usage{
dmaxFratio(x, df, k, log = FALSE)

pmaxFratio(q, df, k, lower.tail = TRUE, log.p = FALSE)

qmaxFratio(p, df, k, lower.tail = TRUE, log.p = FALSE)

rmaxFratio(n, df, k)

smaxFratio(df, k)
}
\arguments{
  \item{x,q}{vector of non-negative quantities}

  \item{p}{vector of probabilities}

  \item{n}{number of values to generate. If n is a vector,
  length(n) values will be generated}

  \item{df}{vector non-negative, integer degrees of
  freedom}

  \item{k}{vector non-negative, integer number of mean
  squares}

  \item{log,log.p}{logical vector; if TRUE, probabilities p
  are given as log(p)}

  \item{lower.tail}{logical vector; if TRUE (default),
  probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X >
  x]}}
}
\value{
The output values conform to the output from other such
functions in \R. \code{dmaxFratio()} gives the density,
\code{pmaxFratio()} the distribution function and
qmaxFratio its inverse. \code{rmaxFratio()} generates
random numbers. \code{smaxFratio()} produces a list
containing parameters corresponding to the arguments --
mean, median, mode, variance, sd, third cental moment,
fourth central moment, Pearson's skewness, skewness, and
kurtosis.
}
\description{
Density, distribution function, quantile function, random
generator and summary function for the maximum F-ratio.
}
\details{
The maximum F-ratio is the ratio of the largest to the
smallest of k independent mean squares, all with the same
df. The usual use is to test for homogeneity of normal
variances.
}
\note{
#' The maximum F-ratio was introduced by Hartley (1950) as
a shortcut test of the homogeneity of variances from normal
data.  Given a set of k mean squares, each based on the
same number of degrees of freedom, df, the test statistic
is the ratio of the largest to the smallest.  Several
tables have been constructed.  The first by David, H.A.
(1952).  Currently the most extensive are those by Nelson
(1987).

It is important to note that tests of this sort are
substantially dependent on the assumption of normality, and
can not be used robustly as can variance ratios in the
analysis of variance.
}
\section{Limitations}{
  The literature contains no information on numerical
  procedures for this distribution, with the result that
  all calculations are slow.

  Finding p from x should give results for almost any
  values of df and k -- of course absolutely enormous
  values will take a while.

  Finding x from p is an iterative calculation dependent on
  a good starting guess.  Such good guesses have been made
  for \eqn{df \le 24}{df <= 24} and \eqn{k \le 160}{df <=
  160}.  NA will be returned if larger values are
  attempted.
}
\examples{
pmaxFratio(4, 10, 10)
pmaxFratio(c(2.3, 4, 8.5), 10, 10)  ## approximately 5\% 50\% and 95\%
qmaxFratio(p=.95,df=c(10,20), k=10)
smaxFratio(10, 10) ## Wait for this, it may take a while
plot(function(x)dmaxFratio(x, 10, 10),1,10)
}
\author{
Bob Wheeler \email{bwheelerg@gmail.com}
}
\references{
Hartley, H.O. (1950) The maximum F-ratio as a short cut
test for heterogeneity of variance. \emph{Biometrika.}
\bold{37.} 308-312.

David, H.A. (1952). Upper 5 and 1\% points of the maximum
F-ratio. \emph{Biometrika.} \bold{38.} 422-424.

Nelson, L.S. (1987). Upper 10\%, 5\% and 1\% points of the
maximum F-ratio, \emph{Jour. Qual. Tech.} \bold{19-3.}
165-167.
}
\keyword{distribution}

