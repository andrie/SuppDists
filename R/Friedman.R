#' Friedman's chi-square
#' 
#' Density, distribution function, quantile function, random generator and summary function for Friedman's chi square.
#' 
#' The Friedman chi-squared is used for nonparametric ANOVA. The data in N rows of an \eqn{N \times r}{N x r} table are ranked separately such that the ranks take the values from 1 to r in the N different rows. The distributions
#' are obtained on the assumption that there is no relationship between the N rows.
#' 
#' Formulae:
#' 
#' Let \eqn{R_j}{Rj} be the sum of ranks for treatment \eqn{j (j=1\dots r)}{j
#' (j=1 \ldots r)}, then the Friedman statistic is
#' 
#' \deqn{ }{x=[12/(N r (r+1))]Sum(1 \ldots r)(Rj^2) - 3 N
#' (r+1)}\deqn{x=\frac{12}{N r (r+1)}\sum_{j=1}^{r}R_j^2 -3N(r+1)}{x=[12/(N r
#' (r+1))]Sum(1 \ldots r)(Rj^2) - 3 N (r+1)}
#' 
#' this is asymptotically equivalent to a \eqn{\chi^2}{chi squared} random
#' variable. One may also calculate the chi squared statistic for the usual
#' analysis of variance which gives
#' 
#' \deqn{ }{F=((N-1) x)/(N (r-1)-x)}\deqn{F=\frac{(N-1)x}{N(r-1)-x}}{F=((N-1)
#' x)/(N (r-1)-x)}
#' 
#' which may be used with the F distribution functions in R for degrees of
#' freedom \eqn{(r-1)} and \eqn{(N-1)(r-1)}.
#' 
#' @aliases dFriedman pFriedman qFriedman rFriedman sFriedman
#' @param x,q vector of non-negative quantities
#' @param p vector of probabilities
#' @param n number of values to generate. If n is a vector, `length(n)` values will be generated
#' @param r vector of number of treatments
#' @param N (N >= 2) vector of number of replications of each treatment
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' 
#' @return The output values conform to the output from other such functions in `Friedman()` gives the density, `pFriedman()` the distribution function and `qFriedman()`its inverse.  `rFriedman()` generates random numbers. `sFriedman()` produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#' @note Exact calculations are made for the following values of the parameters:
#' 
#' \tabular{ll}{ r \tab N \cr 2 \tab 100 \cr 3 \tab 30 \cr 4 \tab 15 \cr 5 \tab
#' 8 } These exact calculations are made using the algorithm of Kendall and
#' Smith (1939).
#' 
#' The incomplete beta, with continuity correction, is used for calculations
#' outside these ranges.  Some appreciation for the accuracy of the
#' approximation may be obtained by comparing the calculated values with exact
#' tables such as Odeh (1977).  Iman and Davenport (1980) have studied the
#' accuracy of various approximations.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references 
#' Kendall, M. and Smith, B.B. (1939). The problem of m rankings.
#' *Ann. Math. Stat.* **10.** 275-287.
#' 
#' Iman, R.L. and Davenport, J.M. (1980). Approximations of the critical region of the Friedman statistic. 
#' *Comm. Stat. Theor. Meth.* **A9(6).** 571-595.
#' 
#' Odeh, R.E. (1977). Extended tables of the distribution of Friedman's S-statistic in the two-way layout. 
#' *Commun. Statist.-Simula. Computa.* **B6(1).** 29-48.
#' @keywords distribution
#' @examples
#' pFriedman(2, r=5, N=10)
#' pFriedman(c(.8, 3.5, 9.3), r=5, N=10) ## approximately 5% 50% and 95%
#' sFriedman(r=5, N=10)
#' plot(function(x) dFriedman(x, r=5, N=10), 0, 10)

#' @export
dFriedman <- function (x, r, N, log = FALSE) {
  M <- max(length(x), length(r), length(N))
  x <- rep(x, length.out = M)
  r <- rep(r, length.out = M)
  N <- rep(N, length.out = M)
  rho <- rep(FALSE, length.out = M)
  value <- .C("dFriedmanR", as.double(x), as.integer(r), as.integer(N), 
              as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}


#' @export
#' @rdname dFriedman
pFriedman <- function (q, r, N, lower.tail = TRUE, log.p = FALSE) {
  M <- max(length(q), length(r), length(N))
  q <- rep(q, length.out = M)
  r <- rep(r, length.out = M)
  N <- rep(N, length.out = M)
  rho <- rep(FALSE, length.out = M)
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
#' @rdname dFriedman
qFriedman <- function (p, r, N, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(r), length(N))
    p <- rep(p, length.out = M)
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    .C("qFriedmanR", as.double(p), as.integer(r), as.integer(N), 
       as.integer(M), as.integer(rho), val = double(M),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dFriedman
rFriedman <- function (n, r, N) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- max(length(r), length(N))
    r <- rep(r, length.out = M)
    N <- rep(N, length.out = M)
    rho <- rep(FALSE, length.out = M)
    .C("rFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
       as.integer(n), as.integer(M), value = double(n),PACKAGE="SuppDists")$value
  }


#' @export
#' @rdname dFriedman
sFriedman <- function (r, N) {
  M <- max(length(r), length(N))
  r <- rep(r, length.out = M)
  N <- rep(N, length.out = M)
  rho <- rep(FALSE, length.out = M)
  value <- .C("sFriedmanR", as.integer(r), as.integer(N), as.integer(rho), 
              as.integer(M), mn = double(M), med = double(M), mod = double(M), 
              var = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
  aList <- list(title = "Friedman's chi-square", r = r, N = N)
  makeStatList(aList, value$mn, value$med, value$var, value$mod, 
               value$third, value$fourth, -1)
}
