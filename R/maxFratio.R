#' The maximum F-ratio distribution.
#' 
#' Density, distribution function, quantile function, random generator and summary function for the maximum F-ratio.
#' 
#' The maximum F-ratio is the ratio of the largest to the smallest of k independent mean squares, all with the same df. The usual use is to test for homogeneity of normal variances.
#' 
#' @aliases maxFratio dmaxFratio pmaxFratio qmaxFratio rmaxFratio smaxFratio
#' @param x,q vector of non-negative quantities
#' @param p vector of probabilities
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param df vector non-negative, integer degrees of freedom
#' @param k vector non-negative, integer number of mean squares
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @return The output values conform to the output from other such functions in \R. \code{dmaxFratio()} gives the density, \code{pmaxFratio()} the distribution function and qmaxFratio its inverse. \code{rmaxFratio()} generates random numbers. \code{smaxFratio()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#' @note
#' #' The maximum F-ratio was introduced by Hartley (1950) as a shortcut test of the homogeneity of variances from normal data.  Given a set of k mean squares, each based on the same number of degrees of freedom, df, the test statistic is the ratio of the largest to the smallest.  Several tables have been constructed.  The first by David, H.A. (1952).  Currently the most extensive are those by Nelson (1987).
#' 
#' It is important to note that tests of this sort are substantially dependent on the assumption of normality, and can not be used robustly as can variance ratios in the analysis of variance.
#' @section Limitations:
#' 
#' The literature contains no information on numerical procedures for this distribution, with the result that all calculations are slow.
#' 
#' Finding p from x should give results for almost any values of df and k -- of course absolutely enormous values will take a while.
#' 
#' Finding x from p is an iterative calculation dependent on a good starting guess.  Such good guesses have been made for \eqn{df \le 24}{df <= 24} and \eqn{k \le 160}{df <= 160}.  NA will be returned if larger values are attempted.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references 
#' Hartley, H.O. (1950) The maximum F-ratio as a short cut test for heterogeneity of variance. 
#' \emph{Biometrika.} \bold{37.} 308-312.
#' 
#' David, H.A. (1952). Upper 5 and 1\% points of the maximum F-ratio.
#' \emph{Biometrika.} \bold{38.} 422-424.
#' 
#' Nelson, L.S. (1987). Upper 10\%, 5\% and 1\% points of the maximum F-ratio,
#' \emph{Jour. Qual. Tech.} \bold{19-3.} 165-167.
#' @keywords distribution
#' @examples
#' 
#' 
#' pmaxFratio(4, 10, 10)
#' pmaxFratio(c(2.3, 4, 8.5), 10, 10)  ## approximately 5% 50% and 95% 
#' qmaxFratio(p=.95,df=c(10,20), k=10)
#' smaxFratio(10, 10) ## Wait for this, it may take a while
#' plot(function(x)dmaxFratio(x, 10, 10),1,10)
#' 
#' @export

dmaxFratio <- function (x, df, k, log = FALSE) {
  if (log == TRUE) 
    p <- exp(p)
  N <- max(length(x), length(df), length(k))
  x <- rep(x, length.out = N)
  df <- rep(df, length.out = N)
  k <- rep(k, length.out = N)
  .C("dmaxFratioR", as.double(x), as.integer(df), as.integer(k), 
     as.integer(N), val = double(N),PACKAGE="SuppDists")$val
}


#' @export
#' @rdname dmaxFratio
pmaxFratio <- function (q, df, k, lower.tail = TRUE, log.p = FALSE) {
    N <- max(length(q), length(df), length(k))
    q <- rep(q, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    if (lower.tail == TRUE) {
      value <- .C("pmaxFratioR", as.double(q), as.integer(df), 
                  as.integer(k), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    else {
      value <- .C("umaxFratioR", as.double(q), as.integer(df), 
                  as.integer(k), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
      value <- log(value)
    value
  }


#' @export
#' @rdname dmaxFratio
qmaxFratio <- function (p, df, k, lower.tail = TRUE, log.p = FALSE) {
    if (lower.tail == FALSE) 
      p <- 1 - p
    if (log.p == TRUE) 
      p <- exp(p)
    N <- max(length(p), length(df), length(k))
    p <- rep(p, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    .C("qmaxFratioR", as.double(p), as.integer(df), as.integer(k), 
       as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  }

#' @export
#' @rdname dmaxFratio
rmaxFratio <- function (n, df, k) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- max(length(df), length(k))
    df <- rep(df, length.out = M)
    k <- rep(k, length.out = M)
    .C("rmaxFratioR", as.integer(df), as.integer(k), as.integer(n), 
       as.integer(M), value = double(n),PACKAGE="SuppDists")$value
  }

#' @export
#' @rdname dmaxFratio
smaxFratio <- function (df, k) {
    N <- max(length(df), length(k))
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    value <- .C("smaxFratioR", as.integer(df), as.integer(k), 
                as.integer(N), mn = double(N), med = double(N), mod = double(N), 
                var = double(N), third = double(N), fourth = double(N),PACKAGE="SuppDists")
    aList <- list(title = "Maximum F ratio", df = df, k = k)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
                 value$third, value$fourth, 2)
  }
