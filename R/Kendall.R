#' The distribution of Kendall's tau.
#' 
#' Density, distribution function, quantile function, random generator and summary function for Kendall's tau.
#' 
#' 
#' There are two categories with N treatments each.  The treatments are ranked for each category, and then sorted according to the ranks for the first category.  This produces a 2 by N array in which the numbers in the first row are increasing from 1 to N.  The array is scanned, and every time two adjacent ranks in the second row are not in order, they are exchanged.  The scanning is repeated until the second row is in increasing order.  Let s denote the number of exchanges, then Kendall's tau is given by
#' 
#' \deqn{\tau=1-\frac{4s}{N(N-1)}}{tau=1-4 s/(N (N-1))}
#' 
#' This too is a product-moment correlation coefficient.  See Kendall (1975), Chapter 2.  Other methods for calculating the statistic are also discussed there.
#' 
#' The calculated values are exact for \eqn{N < 13}, thereafter an Edgeworth expansion is used.
#' 
#' @aliases dKendall pKendall qKendall rKendall sKendall
#' @param x,q vector of non-negative quantities
#' @param p vector of probabilities
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param N vector number of treatments
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @return The output values conform to the output from other such functions in \code{dKendall()} gives the density, \code{pKendall()} the distribution function and \code{qKendall()} its inverse. \code{rKendall()} generates random numbers. \code{sKendall()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references 
#' Kendall, M. (1975). 
#' \emph{Rank Correlation Methods.} Griffin, London.
#' @keywords distribution
#' @examples
#' 
#' pKendall(0, N=10)
#' pKendall(c(-.42, 0.02, .42), N=10) ## approximately 5% 50% and 95% 
#' qKendall(.95, N=c(10, 20))
#' sKendall(N=10)
#' plot(function(x) dKendall(x, N=10), -0.5, 0.5)
#' 
#' 
#' @export
dKendall <- function (x, N, log = FALSE) {
  M <- max(length(x), length(N))
  x <- rep(x, length.out = M)
  N <- rep(N, length.out = M)
  value <- .C("dKendallR", as.integer(N), as.double(x), as.integer(M), 
              val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}


#' @export
#' @rdname dKendall
pKendall <- function (q, N, lower.tail = TRUE, log.p = FALSE) {
  M <- max(length(q), length(N))
  q <- rep(q, length.out = M)
  N <- rep(N, length.out = M)
  if (lower.tail == TRUE) {
    value <- .C("pKendallR", as.integer(N), as.double(q), 
                as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }
  else {
    value <- .C("uKendallR", as.integer(N), as.double(q), 
                as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }
  if (log.p == TRUE) 
    value <- log(value)
  value
}


#' @export
#' @rdname dKendall
qKendall <- function (p, N, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(N))
    p <- rep(p, length.out = M)
    N <- rep(N, length.out = M)
    .C("qKendallR", as.integer(N), as.double(p), as.integer(M), 
       val = double(M),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dKendall
rKendall <- function (n, N) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- length(N)
    .C("rKendallR", as.integer(N), as.integer(n), as.integer(M), 
       val = double(n),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dKendall
sKendall <- function (N) {
    M <- length(N)
    mn <- rep(0, length.out = M)
    med <- rep(0, length.out = M)
    mod <- rep(0, length.out = M)
    third <- rep(0, length.out = M)
    var <- (4 * N + 10)/(9 * N * (N - 1))
    fourth <- .C("fourthKendallR", as.integer(N), as.integer(M), 
                 val = double(M),PACKAGE="SuppDists")$val
    aList <- list(title = "Kendall's Tau", N = N)
    makeStatList(aList, mn, med, var, mod, third, fourth, -1)
  }
