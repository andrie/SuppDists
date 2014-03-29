#' The Pearson product moment correlation coefficient.
#' 
#' Density, distribution function, quantile function, random generator and summary function for the distribution of Pearson's product moment correlation.
#' 
#' 
#' @aliases Pearson dPearson pPearson qPearson rPearson sPearson
#' @param x,q vector of sample correlations
#' @param p vector of probabilities
#' @param rho vector of population correlations
#' @param N vector of numbers of observations, \eqn{(N > 3)}
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are \eqn{P[R <= r]}, otherwise, \eqn{P[R > r]}
#' @return The output values conform to the output from other such functions in \R. \code{dPearson()} gives the density, \code{pPearson()} the distribution function and \code{qPearson()} its inverse. \code{rPearson()} generates random numbers. \code{sPearson()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @keywords distribution
#' @examples
#' 
#' pPearson(0.5, N=10)
#' pPearson(q=0.5, N=10, rho=0.3) 
#' sPearson(N=10)
#' plot(function(x) dPearson(x, N=10, rho=0.7), -1, 1)
#' 
#' @export

dPearson <- function (x, N, rho = 0, log = FALSE) {
  M <- max(length(x), length(rho), length(N))
  x <- rep(x, length.out = M)
  rho <- rep(rho, length.out = M)
  N <- rep(N, length.out = M)
  value <- .C("dcorrR", as.double(x), as.double(rho), as.integer(N), 
              as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}


#' @export
#' @rdname dPearson
pPearson <- function (q, N, rho = 0, lower.tail = TRUE, log.p = FALSE) {
    M <- max(length(q), length(rho), length(N))
    q <- rep(q, length.out = M)
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    if (lower.tail == TRUE) {
      value <- .C("pcorrR", as.double(q), as.double(rho), as.integer(N), 
                  as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    else {
      value <- .C("ucorrR", as.double(q), as.double(rho), as.integer(N), 
                  as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
      value <- log(value)
    value
  }


#' @export
#' @rdname dPearson
qPearson <- function (p, N, rho = 0, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(rho), length(N))
    p <- rep(p, length.out = M)
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    .C("qcorrR", as.double(p), as.double(rho), as.integer(N), 
       as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dPearson
rPearson <- function (n, N, rho = 0) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- max(length(rho), length(N))
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    .C("rcorrR", as.double(rho), as.integer(N), as.integer(n), 
       as.integer(M), val = double(n),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dPearson
sPearson <- function (N, rho = 0) {
    M <- max(length(rho), length(N))
    rho <- rep(rho, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("scorrR", as.double(rho), as.integer(N), as.integer(M), 
                mn = double(M), med = double(M), mod = double(M), var = double(M), 
                third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Correlation coefficient", rho = rho, 
                  N = N)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
                 value$third, value$fourth, -1)
  }
