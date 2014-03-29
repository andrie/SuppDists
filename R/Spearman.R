#' Spearman's rho.
#' 
#' Density, distribution function, quantile function, random generator and summary function for Spearman's rho.
#' 
#' 
#' Spearman's rho is the rank correlation coefficient between r pairs of items. It ranges from -1 to 1. Denote by d, the sum of squares of the differences between the matched ranks, then x is given by:
#' 
#' \deqn{1-\frac{6d}{r(r^2-1)}}{1-6 d /(r(r^2-1))}
#' 
#' This is, in fact, the product-moment correlation coefficient of rank differences.  See Kendall (1975), Chapter 2. It is identical to Friedman's chi-squared for two treatments scaled to the -1, 1 range -- if X is the Friedman statistic, then \eqn{\rho=frac{X}{r-1)-1}}{rho = X/(r-1) -1}.
#' 
#' Exact calculations are made for \eqn{r \le 100}{r <= 100}
#' 
#' These exact calculations are made using the algorithm of Kendall and Smith (1939).
#' 
#' The incomplete beta, with continuity correction, is used for calculations outside this range.
#' 
#' @aliases dSpearman pSpearman qSpearman rSpearman sSpearman
#' @param x,q vector of non-negative quantities
#' @param p vector of probabilities
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param r (r >= 3) vector of number of observations
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @return The output values conform to the output from other such functions in \R. \code{dSpearman()} gives the density, \code{pSpearman()} the distribution function and \code{qSpearman()} its inverse. \code{rSpearman()} generates random numbers. \code{sSpearman()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references
#' 
#' Kendall, M. (1975). \emph{Rank Correlation Methods.} Griffin, London.
#' 
#' Kendall, M. and Smith, B.B. (1939). The problem of m rankings. 
#' \emph{Ann. Math. Stat.} \bold{10.} 275-287.
#' @keywords distribution
#' @examples
#' 
#' pSpearman(.95, 10)
#' pSpearman(c(-0.55, 0, 0.55), 10) ## approximately 5% 50% and 95% 
#' sSpearman(10)
#' plot(function(x) dSpearman(x, 10), -.9, .9)
#' 
#' @export
dSpearman <- function (x, r, log = FALSE) {
  M <- max(length(x), length(r))
  x <- rep(x, length.out = M)
  r <- rep(r, length.out = M)
  N <- rep(2, length.out = M)
  rho <- rep(TRUE, length.out = M)
  value <- .C("dFriedmanR", as.double(x), as.integer(r), as.integer(N), 
              as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}

#' @export
#' @rdname dSpearman
pSpearman <- function (q, r, lower.tail = TRUE, log.p = FALSE) {
    M <- max(length(q), length(r))
    q <- rep(q, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    if (lower.tail == TRUE) {
      value <- .C("pFriedmanR", as.double(q), as.integer(r), 
                  as.integer(N), as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
      value <- .C("uFriedmanR", as.double(q), as.integer(r), 
                  as.integer(N), as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
      value <- log(value)
    value
  }


#' @export
#' @rdname dSpearman
qSpearman <- function (p, r, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(r))
    p <- rep(p, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    .C("qFriedmanR", as.double(p), as.integer(r), as.integer(N), 
       as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dSpearman
rSpearman <- function (n, r) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- length(r)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    .C("rFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
       as.integer(n), as.integer(M), value = double(n),PACKAGE="SuppDists")$value
  }

#' @export
#' @rdname dSpearman
sSpearman <- function (r) {
    M <- length(r)
    r <- rep(r, length.out = M)
    N <- rep(2, length.out = M)
    rho <- rep(TRUE, length.out = M)
    value <- .C("sFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
                as.integer(M), mn = double(M), med = double(M), mod = double(M), 
                var = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Spearman's rho", r = r)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
                 value$third, value$fourth, -1)
  }
