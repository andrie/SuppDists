\name{ghyper}
\alias{dghyper}
\alias{pghyper}
\alias{qghyper}
\alias{rghyper}
\alias{sghyper}
\alias{tghyper}
\alias{Generalized hypergeometric}
\alias{Negative hypergeometric}
\alias{Inverse hypergeometric}
\alias{Hypergeometric waiting time}
\alias{Beta-binomial}
\alias{Beta-negative-binomial}
\alias{Beta-Pascal}
\alias{Generalized Waring} 


\title{Generalized hypergeometric distributions}
\description{
  Density, distribution function, quantile function, random generator and summary function for generalized hypergeometric distributions.
}
\usage{
dghyper(x, a, k, N, log=FALSE)
pghyper(q, a, k, N, lower.tail=TRUE, log.p=FALSE)
qghyper(p, a, k, N, lower.tail=TRUE, log.p=FALSE)
rghyper(n, a, k, N)
sghyper(a, k, N)
tghyper(a, k, N) ## scalar arguments only
}
\arguments{
	\item{x,q,n}{vector of non-negative \bold{integer} quantities}
	\item{p}{vector of probabilities}
	\item{a}{vector of real values giving the first column total}
	\item{k}{vector  of real values giving the first row total}
	\item{N}{vector of real values giving the grand total}
	\item{log, log.p}{logical vector; if TRUE, probabilities p are given as log(p)}
	\item{lower.tail}{logical vector; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}}
}
\value{
  The output values conform to the output from other such functions in R. \code{dghyper()} gives the density, \code{pghyper()} the distribution function and \code{qghyper()} its inverse. \code{rghyper()} generates random numbers. \code{sghyper()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third central moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
\cr
	The function \code{tghyper()} returns the hypergeometric type and the range of values for x.	
}


\details{

The basic representation is in terms of a two-way table:

\tabular{ccc}{
x \tab  k-x \tab  k\cr
a-x \tab  b-k+x \tab  N-k\cr
a \tab     b \tab  N}

and the associated hypergeometric probability \eqn{P(x)=C_x^a C_{k-x}^b / C_k^N}{P(x)=choose(a, x) choose(b, k-x) / choose(N, k)}.

The table is constrained so that rows and columns add to the margins.  In all cases x is an integer or zero, but meaningful probability distributions occur when the other parameters are real. Johnson, Kotz and Kemp (1992) give a general discussion.

Kemp and Kemp (1956) classify the possible probability distributions that can occur when real values are allowed, into eight types. The \emph{classic} hypergeometric with integer values forms a ninth type.  Five of the eight types correspond to known distributions used in various contexts.  Three of the eight types, appear to have no practical applications, but for completeness they have been implemented.  

The Kemp and Kemp types are defined in terms of the ranges of the a, k, and N parameters and are given in \code{\link{ghyper.types}}. The function \code{tghyper()} will give details for specific values of a, k, and N.
 


These distributions apply to many important problems, which has lead to a variety of names:

The Kemp and Kemp types IIA and IIIA are known as:
\itemize{
\item Negative hypergeometric
\item	Inverse hypergeometric
\item	Hypergeometric waiting time
\item	Beta-binomial
}

The advantages of the conditional argument are considerable. Consider a few examples:

\enumerate{
\item Future event: Consider two events which have occurred u and v times respectively.  The distribution function of x occurrences of the first event in a sample of k new trials is calculated.  Here a = -u-1, and N = -u-v-2. 

Example: Suppose Toronto has won 3 games and Atlanta 1 in the World Series. What is the probability that Toronto will win the series by taking 2 or more of the remaining 3 games?

\item Exceedance: Consider two samples of size m and k, then the distribution function of x, the number of elements out of k which exceed the r th largest element in the size m sample is calculated.  Here a = -r, and N = -m-1. 

Example: Suppose that only once in the last century has the high-water mark at the St. Joe bridge exceeded 12 feet, what is the probability that it will not do so in the next ten years?

\item Waiting time: Consider an urn with T balls, m of which are white, and that drawing without replacement is continued until w white balls are obtained, then the distribution function of x, the number of balls in excess of w that must be drawn is desired.  Here a = -w , N = -m-1, and k = T - m. 

Example: Suppose a lot of 100 contains 5 defectives. What is the mean number of items that must be inspected before a defective item is found?

\item Mixture: Suppose x has a binomial distribution with parameter p, and number of trials k.  Suppose that p is not fixed, but itself distributed like a beta variable with parameters A and B, then the distribution of x is calculated with a = -A and N = -A -B.
}
 



Names for Kemp and Kemp type IV are:
\itemize{
\item	Beta-negative-binomial
\item	Beta-Pascal
\item	Generalized Waring 
}

One application is accidents:

Suppose accidents follow a Poisson distribution with mean L, and suppose L varies with individuals according to accident proneness, m.  In particular, suppose L follows a gamma distribution with parameter r and scale factor m , and that the scale factor n itself follows a beta distribution with parameters A and B, then the distribution of accidents, x, is beta-negative-binomial with a = -B, k = -r , and N = A -1.  See Xekalki (1983) for a discussion of this as well as a discussion of accident models for proneness, contagion and spells.



}

\note{The parameters of these functions differ from those of the
  hypergeometric functions of R. To translate between the two use the
  following as a model: phyper(x,m,n,k) = pghyper(x,k,m,m+n).
}


\references{

Johnson, N.L., Kotz, S. and Kemp, A. (1992) \emph{Univariate discrete distributions}. Wiley, N.Y.
 
Kemp, C.D., and Kemp, A.W. (1956). Generalized hypergeometric distributions. \emph{Jour. Roy. Statist. Soc. B.} \bold{18.} 202-211.



Xekalaki, E. (1983). The univariate generalized Waring distribution in relation to accident theory: proneness, spells or contagion. \emph{Biometrics.} \bold{39.}  887-895.
}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
}
\examples{

tghyper(a=4, k=4, N=10) 		## classic
tghyper(a=4.1, k=5, N=10) 		## type IA(i) Real classic
tghyper(a=5, k=4.1, N=10) 		## type IA(ii) Real classic
tghyper(a=4.2, k=4.6, N=12.2) 		## type IB
tghyper(a=-5.1, k=10, N=-7) 		## type IIA
tghyper(a=-0.5, k=5.9, N=-0.7) 		## type IIB
tghyper(a=10, k=-5.1, N=-7) 		## type IIIA Negative hypergeometric
tghyper(a=5.9, k=-0.5, N=-0.7) 		## type IIIB
tghyper(a=-1, k=-1, N=5) 		## type IV Generalized Waring

sghyper(a=-1, k=-1, N=5)
plot(function(x)dghyper(x,a=-1,k=-1,N=5),0,5)

#Fisher's exact test: contingency table with rows (1,3),(3,1) 
pghyper(1,4,4,8)
pghyper(3,4,4,8,lower.tail=FALSE)



#Beta-binomial applications:

#Application examples:
tghyper(-4,3,-6)
pghyper(2,-4,3,-6,lower=FALSE)
pghyper(0,-2,10,-101)
sghyper(-1,95,-6)$Mean+1



}
\keyword{distribution}