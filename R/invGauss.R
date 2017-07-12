#' The inverse Gaussian and Wald distributions.
#' 
#' Density, distribution function, quantile function, random generator and summary function for the inverse Gaussian and Wald distributions.
#' 
#' Probability functions:
#' 
#' \deqn{f(x,\nu,\lambda)=\sqrt{\frac{\lambda}{2\pi x^3}}\exp\left[-\lambda\frac{(x-\nu)^2}{2\nu^2 x}\right]
#' \mbox{-- the density}}{f(x,nu,lambda)=sqrt[lambda/(2 pi x^3)]exp[-lambda(x-nu)^2/(2 x nu^2)] -- the density }
#' 
#' \deqn{F(x,\nu,\lambda)=\Phi\left[\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\nu}-1\right)\right]+e^{2\lambda/\nu}\Phi\left[\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\nu}+1\right)\right]
#' \mbox{-- the distribution function}}{F(x,nu,lambda)=Phi[sqrt(lambda/x)((x/nu)-1)]+exp(2*lambda/nu)*Phi[sqrt(lambda/x)((x/nu)+1)] -- the distribution function}
#' 
#' where \eqn{\Phi[]}{Phi[]} is the standard normal distribution function.
#' 
#' The calculations are accurate to at least seven significant figures over an extended range - much larger than that of any existing tables. We have tested them for \eqn{\lambda / \nu = 10e-20}{(lambda/nu)=10e-20}, and \eqn{\lambda/\nu=10e4}{lambda / nu = 10e4}. Accessible tables are those of Wassan and Roy (1969), which unfortunately, are sometimes good to only two significant digits. Much better tables are available in an expensive CRC Handbook (1989), which are accurate to at least 7 significant digits for \eqn{\lambda/\nu \ge 0.02}{lambda / nu >= 0.02} to \eqn{\lambda/\nu \le 4000}{lambda / nu <= 4000}.
#' 
#' These are first passage time distributions of Brownian motion with positive drift. See Whitmore and Seshadri (1987) for a heuristic derivation. The Wald (1947) form represents the average sample number in sequential analysis. The distribution has a non-monotonic failure rate, and is of considerable interest in lifetime studies: Ckhhikara and Folks (1977). A general reference is Seshadri (1993).
#' 
#' This is an extremely difficult distribution to treat numerically, and it would not have been possible without some extraordinary contributions. An elegant derivation of the distribution function is to be found in Shuster (1968). The first such derivation seems to be that of Zigangirov (1962), which because of its inaccessibility, the author has not read. The method of generating random numbers is due to Michael, Schucany, and Haas (1976). The approximation of Whitmore and Yalovsky (1978) makes it possible to find starting values for inverting the distribution. All three papers are short, elegant, and non- trivial.
#' 
#' @aliases inverse-Gaussian Wald-distribution dinvGauss pinvGauss qinvGauss rinvGauss sinvGauss
#' @param x,q vector of non-negative quantities
#' @param p vector of probabilities
#' @param n vector of numbers of observations
#' @param nu vector real and non-negative parameter -- the Wald distribution results when nu=1
#' @param lambda vector real and non-negative parameter
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @return The output values conform to the output from other such functions in R. `dinvGauss()` gives the density, `pinvGauss()` the distribution function and `qinvGauss()` its inverse. `rinvGauss()` generates random numbers. `sinvGauss()` produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references 
#' Ckhhikara, R.S. and Folks, J.L. (1977) The inverse Gaussian distribution as a lifetime model. 
#' \emph{Technometrics.} \bold{19-4.} 461-468.
#' 
#' CRC Handbook. (1989). 
#' \emph{Percentile points of the inverse Gaussian distribution.} J.A. Koziol (ed.) Boca Raton, FL.
#' 
#' Michael, J.R., Schucany, W.R. and Haas, R.W. (1976). Generating random variates using transformations with multiple roots. 
#' \emph{American Statistician.} \bold{30-2.} 88-90.
#' 
#' Seshadri, V. (1993). 
#' \emph{The inverse Gaussian distribution.} Clarendon, Oxford
#' 
#' Shuster, J. (1968). On the inverse Gaussian distribution function. Jour.
#' \emph{Am. Stat. Assoc.} \bold{63.} 1514-1516.
#' 
#' Wasan, M.T. and Roy, L.K. (1969). Tables of inverse Gaussian percentage points. 
#' \emph{Technometrics.} \bold{11-3.} 591-604.
#' 
#' Wald, A. (1947). 
#' \emph{Sequential analysis.} Wiley, NY
#' 
#' Whitmore, G.A. and Seshadri, V. (1987). A heuristic derivation of the inverse Gaussian distribution. 
#' \emph{American Statistician.} \bold{41-4.} 280-281.
#' 
#' Whitmore, G.A. and Yalovsky, M. (1978). A normalizing logarithmic transformation for inverse Gaussian random variables. 
#' \emph{Technometrics.} \bold{20-2.} 207-208.
#' 
#' Zigangirov, K.S. (1962). Expression for the Wald distribution in terms of normal distribution. 
#' \emph{Radiotech.Electron.} \bold{7.} 164-166.
#' @keywords distribution
#' @examples
#' 
#' pinvGauss(1, 1, 16)
#' pinvGauss(c(.65, 1, 1.45), 1, 16) ## approximately 5% 50% and 95%
#' pars <- sinvGauss(1, 16)
#' plot(function(x)dinvGauss(x,1, 16), pars$Mean - 3*pars$SD, pars$Mean + 3*pars$SD)
#' 

#' @export
dinvGauss <- function (x, nu, lambda, log = FALSE) {
  N <- max(length(x), length(nu), length(lambda))
  x <- rep(x, length.out = N)
  nu <- rep(nu, length.out = N)
  lambda <- rep(lambda, length.out = N)
  value <- .C("dinvGaussR", as.double(x), as.double(nu), as.double(lambda), 
              as.integer(N), lambda = double(N),PACKAGE="SuppDists")$lambda
  if (log == TRUE) 
    value <- log(value)
  value
}

#' @export
#' @rdname dinvGauss
pinvGauss <- function (q, nu, lambda, lower.tail = TRUE, log.p = FALSE) {
  N <- max(length(q), length(nu), length(lambda))
  q <- rep(q, length.out = N)
  nu <- rep(nu, length.out = N)
  lambda <- rep(lambda, length.out = N)
  if (lower.tail == TRUE) {
    value <- .C("pinvGaussR", as.double(q), as.double(nu), 
                as.double(lambda), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  }
  else {
    value <- .C("uinvGaussR", as.double(q), as.double(nu), 
                as.double(lambda), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  }
  if (log.p == TRUE) 
    value <- log(value)
  value
}


#' @export
#' @rdname dinvGauss
qinvGauss <- function (p, nu, lambda, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    N <- max(length(p), length(nu), length(lambda))
    p <- rep(p, length.out = N)
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    .C("qinvGaussR", as.double(p), as.double(nu), as.double(lambda), 
       as.integer(N), value = double(N),PACKAGE="SuppDists")$value
  }

#' @export
#' @rdname dinvGauss
rinvGauss <- function (n, nu, lambda) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    N <- max(length(nu), length(lambda))
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    .C("rinvGaussR", as.double(nu), as.double(lambda), as.integer(n), 
       as.integer(N), value = double(n),PACKAGE="SuppDists")$value
  }


#' @export
#' @rdname dinvGauss
sinvGauss <- function (nu, lambda) {
    N <- max(length(nu), length(lambda))
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    med <- qinvGauss(0.5, nu, lambda)
    nu[nu<=0]<-NA
    lambda[lambda<=0]<-NA
    factor <- (nu^2)/lambda
    var <- nu * factor
    k3 <- 3 * var * factor
    k4 <- 5 * k3 * factor
    mod <- -1.5 * factor + nu * sqrt(1 + 2.25 * (nu/lambda)^2)
    third <- k3
    fourth <- k4 + 3 * var^2
    aList <- list(title = "Inverse Gaussian", nu = nu, lambda = lambda)
    makeStatList(aList, nu, med, var, mod, third, fourth, -1)
  }
