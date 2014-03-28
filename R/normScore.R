#' Normal Scores distribution.
#' 
#' Density, distribution function, quantile function, random generator and summary function for the normal scores test. A function to calculate expected values of normal order statistics is included.
#' 
#' This is the Kruskal-Wallis statistic with ranks replaced by the expected values of normal order statistics. There are c treatments with sample sizes \eqn{n_j, j=1 \dots c}{nj, j=1 \ldots c}. The total sample size is \eqn{N=\sum_1^c n_j}{N=Sum nj}. The distribution depends on c, N, and U, where \eqn{U=\sum_1^c (1/n_j)}{U=Sum (1/nj)}.
#' 
#' Let \eqn{e_N(k)}{eN(k)} be the expected value of the \eqn{k_{th}}{kth} smallest observation in a sample of N independent normal variates. Rank all observations together, and let \eqn{R_{ij}}{Rij} denote the rank of observation \eqn{X_{ij}}{Xij}, \eqn{i=1 \dots n_j}{i=1 \ldots nj} for treatment \eqn{j=1 \dots c}{j=1 \ldots c}, then the normal scores test statistic is
#' 
#' \deqn{x=(N-1)\frac{1}{\sum_{k=1}^{N} e_N(k)^2} \sum_{j=1}^{c}\frac{S_j^2}{n_j}}{x=(N-1)[1/Sum(1 \ldots N)(eN(k)^2)]Sum(1 \ldots c)[(Sj^2)/nj]}
#' 
#' where \eqn{S_j=\sum_{i=1}^{n_j}(e_N(R_{ij}))}{Sj=Sum(1 \ldots nj)(eN(Rij))}.
#' 
#' See Lu and Smith (1979) for a thorough discussion and some exact tables for small r and n.  The calculations made here use an incomplete beta approximation -- the same one used for Kruskal-Wallis, only differing in the calculation of the variance of the statistic.
#' 
#' The expected values of the normal order statistics use a modification of M.Maechler's C version of the Fortran algorithm given by Royston (1982). Spot checking the values against Harter (1969) confirms the accuracy to 4 decimal places as claimed by Royston.
#' 
#' @aliases NormScore dNormScore pNormScore qNormScore rNormScore sNormScore normOrder
#' @param x,q vector of non-negative quantities
#' @param p vector of probabilities
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param c vector number of treatments
#' @param N vector total number of observations
#' @param U vector sum of reciprocals of the number of the c sample sizes
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @return The output values conform to the output from other such functions in \R. \code{dNormScore()} gives the density, \code{pNormScore()} the distribution function and \code{qNormScore()} its inverse. \code{rNormScore()} generates random numbers. \code{sNormScore()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis. \code{normOrder()} gives the expected values of the normal order statistics for a sample of size N.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references
#' 
#' Harter, H.L. (1969). 
#' \emph{Order statistics and their use in testing and estimation, volume 2.} U.S. Supp. of Doc.
#' 
#' Lu, H.T. and Smith, P.J. (1979) Distribution of normal scores statistic for nonparametric one-way analysis of variance. 
#' \emph{Jour. Am Stat. Assoc.} \bold{74.} 715-722.
#' 
#' Royston, J.P. (1982). Expected normal order statistics (exact and approximate) AS 177. 
#' \emph{Applied Statistics.} \bold{31.} 161-165.
#' @keywords distribution
#' @examples
#' 
#' #Assuming three treatments, each with a sample size of 5
#' 
#' pNormScore(2, 3, 15, 0.6)
#' pNormScore(c(0.11,1.5,5.6), 3, 15, 0.6) ## approximately 5% 50% and 95%
#' sNormScore(3, 15, 0.6)
#' plot(function(x)dNormScore(x,c=5, N=15, U=0.6),0,5)
#' 
#' @export
dNormScore <- function (x, c, N, U, log = FALSE) {
  M <- max(length(x), length(c), length(N), length(U))
  x <- rep(x, length.out = M)
  c <- rep(c, length.out = M)
  n <- rep(N, length.out = M)
  U <- rep(U, length.out = M)
  Ns <- rep(TRUE, length.out = M)
  value <- .C("dKruskalWallisR", as.double(x), as.integer(c), 
              as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
              val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}


#' @export
#' @rdname dNormScore
pNormScore <- function (q, c, N, U, lower.tail = TRUE, log.p = FALSE) {
    M <- max(length(q), length(c), length(N), length(U))
    q <- rep(q, length.out = M)
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    if (lower.tail == TRUE) {
      value <- .C("pKruskalWallisR", as.double(q), as.integer(c), 
                  as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
                  val = double(M),PACKAGE="SuppDists")$val
    }
    else {
      value <- .C("uKruskalWallisR", as.double(q), as.integer(c), 
                  as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
                  val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
      value <- log(value)
    value
  }


#' @export
#' @rdname dNormScore
qNormScore <- function (p, c, N, U, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(c), length(N), length(U))
    p <- rep(p, length.out = M)
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    .C("qKruskalWallisR", as.double(p), as.integer(c), as.integer(N), 
       as.double(U), as.integer(Ns), as.integer(M), val = double(M),
       PACKAGE="SuppDists")$val
  }


#' @export
#' @rdname dNormScore
rNormScore <- function (n, c, N, U) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    .C("rKruskalWallisR", randArray = double(n), as.integer(n), 
       as.integer(M), as.integer(c), as.integer(N), as.double(U), 
       as.integer(Ns),PACKAGE="SuppDists" )$randArray
  }


#' @export
#' @rdname dNormScore
sNormScore <- function (c, N, U) {
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(TRUE, length.out = M)
    value <- .C("sKruskalWallisR", as.integer(c), as.integer(n), 
                as.double(U), as.integer(Ns), as.integer(M), var = double(M), 
                mod = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    mn <- (c - 1)
    aList <- list(title = "Normal Scores", c = c, N = n, U = U)
    median <- qNormScore(0.5, c, n, U)
    makeStatList(aList, mn, median, value$var, value$mod, value$third, 
                 value$fourth, -1)
  }
