#' The Johnson distributions.
#' 
#' Density, distribution function, quantile function, random generator and summary function for the Johnson distributions.
#' 
#' The Johnson system (Johnson 1949) is a very flexible system for describing statistical distributions. It is defined by 
#' \deqn{z=\gamma+\delta \log{f(u)}, u=(x-\xi)/\lambda}{z=gamma+delta log(f(u)), with u=(x-xi)/lambda}
#' 
#' and where \eqn{f( )} has four possible forms:
#' 
#' \tabular{ll}{ SL:\tab \eqn{f(u)=u} the log normal \cr SU:\tab
#' \eqn{f(u)=u+\sqrt{1+u^2}}{f(u)=u+sqrt(1+u^2)} an unbounded distribution\cr
#' SB:\tab \eqn{f(u)=u/(1-u)} a bounded distribution\cr SN:\tab \eqn{\exp(u)}
#' the normal }
#' 
#' Estimation of the Johnson parameters may be done from quantiles. The procedure of Wheeler (1980) is used.
#' 
#' They may also be estimated from the moments.  Applied Statistics algorithm 99, due to Hill, Hill, and Holder (1976) has been translated into C for this implementation.
#' 
#' @aliases dJohnson pJohnson qJohnson rJohnson sJohnson moments JohnsonFit
#' @param x,q vector of quantities
#' @param t observation vector, t=x, or moment vector, t=[mean,m2,m3,m4]
#' @param p vector of probabilities
#' @param n vector of numbers of observations
#' @param parms list or list of lists each containing output of JohnsonFit()
#' @param moment character scalar specifying t: "quant" (default), or "use," or "find"
#' @param log,log.p logical vector; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical vector; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @return The output values conform to the output from other such functions in \code{dJohnson()} gives the density, \code{pJohnson()} the distribution function and \code{qJohnson()} its inverse. \code{rJohnson()} generates random numbers. \code{sJohnson()} produces a list containing parameters corresponding to the arguments -- mean, median, mode, variance, sd, third cental moment, fourth central moment, Pearson's skewness, skewness, and kurtosis.
#'
#' \code{moments()} calculates the moment statistics of x as a vector with elements (mu, sigma, skew, kurt), where mu is the mean of x, sigma the SD of x with divisor \code{length(x)}, skew is the skewness and kurt the kurtosis.
#' 
#' \code{JohnsonFit()} outputs a list containing the Johnson parameters (gamma, delta, xi, lambda, type), where type is one of the Johnson types: "SN", "SL", "SB", or "SU". \code{JohnsonFit()} does this using 5 order statistics when moment="quant", when moment="find" it does this by using the first four moments of t calculated by the function \code{moments()}, when moment="use" it assumes that the vector t is [mean,m2,m3,m4], where mi is the ith moment about the mean.
#' Fitting by moments is difficult numerically and often \code{JohnsonFit()} will report an error.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references
#' 
#' Hill, I.D., Hill, R., and Holder, R.L. (1976). Fitting Johnson curves by moments. 
#' \emph{Applied Statistics.} AS99.
#' 
#' Johnson, N.L. (1949). Systems of frequency curves generated by methods of translation. 
#' \emph{Biometrika,} \bold{36.} 149-176.
#' 
#' Wheeler, R.E. (1980). Quantile estimators of Johnson curve parameters.
#' \emph{Biometrika.} \bold{67-3} 725-728
#' @keywords distribution
#' @examples
#' 
#' 
#' xx <- rnorm(500)
#' parms <- JohnsonFit(xx)
#' sJohnson(parms)
#' plot(function(xx) dJohnson(xx, parms), -2, 2)
#' pJohnson(1, parms)
#' parms2 <- JohnsonFit(rexp(50))
#' qJohnson(p=0.5, list(parms, parms2))
#' 
#' ## JohnsonFit with moment="find" and moment="use" is not always possible,
#' ## and even when possible, may produce odd results.
#' ## parms <- JohnsonFit(x,moment="find")
#' 
#' parms <- JohnsonFit(c(0, 1, -.5, 4), moment="use")
#' 
#' sJohnson(parms) 
#' 
#' # Fit illustration
#' data(cars)
#' xx <- cars$speed
#' parms <- JohnsonFit(xx)
#' hist(xx, freq=FALSE)
#' plot(function(x)dJohnson(x, parms), 0, 25, add=TRUE)
#' 
#' 
#' @export
dJohnson <- function (x, parms, log = FALSE) {
  tfun <- function(x) if (x == "SN") 
    1
  else if (x == "SL") 
    2
  else if (x == "SU") 
    3
  else 4
  vecFromList <- function(item, aList) {
    if (!is.list(aList[[1]])) 
      return(aList[[item]])
    else {
      tVec <- vector(length = 0)
      for (i in 1:length(aList)) {
        tVec <- append(tVec, (aList[[i]])[[item]])
      }
    }
    tVec
  }
  gamma <- vecFromList(1, parms)
  delta <- vecFromList(2, parms)
  xi <- vecFromList(3, parms)
  lambda <- vecFromList(4, parms)
  type <- vecFromList(5, parms)
  type <- sapply(type, tfun)
  N <- max(length(gamma), length(x))
  x <- rep(x, length.out = N)
  gamma <- rep(gamma, length.out = N)
  delta <- rep(delta, length.out = N)
  xi <- rep(xi, length.out = N)
  lambda <- rep(lambda, length.out = N)
  type <- rep(type, length.out = N)
  value <- .C("dJohnsonR", as.double(x), as.double(gamma), 
              as.double(delta), as.double(xi), as.double(lambda), as.integer(type), 
              as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}



#' @export
#' @importFrom stats quantile var
#' @rdname dJohnson
JohnsonFit <-   function (t, moment = "quant") {
  firstChar=substring(moment,1,1)
  if (firstChar=="f") {
    mom <- moments(t)
    mu <- mom[[1]]
    sigma <- mom[[2]]
    skew <- mom[[3]]
    kurt <- mom[[4]]
    value <- .C("JohnsonMomentFitR", as.double(mu), as.double(sigma), 
                as.double(skew), as.double(kurt), gamma = double(1), 
                delta = double(1), xi = double(1), lambda = double(1), 
                type = integer(1),PACKAGE="SuppDists")
  }
  else if (firstChar=="u") {
    mu<-t[1]
    sigma<-sqrt(t[2])
    skew<-t[3]/sigma^3
    kurt<-(t[4]/t[2]^2)-3
    value <- .C("JohnsonMomentFitR", as.double(mu), as.double(sigma), 
                as.double(skew), as.double(kurt), gamma = double(1), 
                delta = double(1), xi = double(1), lambda = double(1), 
                type = integer(1),PACKAGE="SuppDists")
  }
  else if (firstChar=="q") {
    input <- quantile(t, probs = c(0.05, 0.206, 0.5, 0.794, 
                                   0.95), names = FALSE)
    x5 <- input[[1]]
    x20.6 <- input[[2]]
    x50 <- input[[3]]
    x79.4 <- input[[4]]
    x95 <- input[[5]]
    value <- .C("JohnsonFitR", as.double(x95), as.double(x79.4), 
                as.double(x50), as.double(x20.6), as.double(x5), 
                gamma = double(1), delta = double(1), xi = double(1), 
                lambda = double(1), type = integer(1),PACKAGE="SuppDists")
  }
  else return(NA)
  types <- c("SN", "SL", "SU", "SB")
  list(gamma = value$gamma, delta = value$delta, xi = value$xi, 
       lambda = value$lambda, type = types[value$type])
}


#' @export
#' @rdname dJohnson
pJohnson <- function (q, parms, lower.tail = TRUE, log.p = FALSE) {
  tfun <- function(x) if (x == "SN") 
    1
  else if (x == "SL") 
    2
  else if (x == "SU") 
    3
  else 4
  vecFromList <- function(item, aList) {
    if (!is.list(aList[[1]])) 
      return(aList[[item]])
    else {
      tVec <- vector(length = 0)
      for (i in 1:length(aList)) {
        tVec <- append(tVec, (aList[[i]])[[item]])
      }
    }
    tVec
  }
  gamma <- vecFromList(1, parms)
  delta <- vecFromList(2, parms)
  xi <- vecFromList(3, parms)
  lambda <- vecFromList(4, parms)
  type <- vecFromList(5, parms)
  type <- sapply(type, tfun)
  N <- max(length(gamma), length(q))
  q <- rep(q, length.out = N)
  gamma <- rep(gamma, length.out = N)
  delta <- rep(delta, length.out = N)
  xi <- rep(xi, length.out = N)
  lambda <- rep(lambda, length.out = N)
  type <- rep(type, length.out = N)
  if (lower.tail == TRUE) {
    value <- .C("pJohnsonR", as.double(q), as.double(gamma), 
                as.double(delta), as.double(xi), as.double(lambda), 
                as.integer(type), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  }
  else {
    value <- .C("uJohnsonR", as.double(q), as.double(gamma), 
                as.double(delta), as.double(xi), as.double(lambda), 
                as.integer(type), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  }
  if (log.p == TRUE) 
    value <- log(value)
  value
}

#' @export
#' @rdname dJohnson
qJohnson <- function (p, parms, lower.tail = TRUE, log.p = FALSE) {
    tfun <- function(x) if (x == "SN") 
      1
    else if (x == "SL") 
      2
    else if (x == "SU") 
      3
    else 4
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    vecFromList <- function(item, aList) {
      if (!is.list(aList[[1]])) 
        return(aList[[item]])
      else {
        tVec <- vector(length = 0)
        for (i in 1:length(aList)) {
          tVec <- append(tVec, (aList[[i]])[[item]])
        }
      }
      tVec
    }
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    N <- max(length(gamma), length(p))
    p <- rep(p, length.out = N)
    gamma <- rep(gamma, length.out = N)
    delta <- rep(delta, length.out = N)
    xi <- rep(xi, length.out = N)
    lambda <- rep(lambda, length.out = N)
    type <- rep(type, length.out = N)
    .C("qJohnsonR", as.double(p), as.double(gamma), as.double(delta), 
       as.double(xi), as.double(lambda), as.integer(type), as.integer(N), 
       val = double(N),PACKAGE="SuppDists")$val
  }


#' @export
#' @rdname dJohnson
rJohnson <- function (n, parms) {
    tfun <- function(x) if (x == "SN") 
      1
    else if (x == "SL") 
      2
    else if (x == "SU") 
      3
    else 4
    vecFromList <- function(item, aList) {
      if (!is.list(aList[[1]])) 
        return(aList[[item]])
      else {
        tVec <- vector(length = 0)
        for (i in 1:length(aList)) {
          tVec <- append(tVec, (aList[[i]])[[item]])
        }
      }
      tVec
    }
    n <- if (length(n) > 1) 
      length(n)
    else n
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    M <- length(gamma)
    .C("rJohnsonR", as.double(gamma), as.double(delta), as.double(xi), 
       as.double(lambda), as.integer(type), as.integer(n), as.integer(M), 
       val = double(n),PACKAGE="SuppDists")$val
  }


#' @export
#' @rdname dJohnson
sJohnson <- function (parms) {
    tfun <- function(x) if (x == "SN") 
      1
    else if (x == "SL") 
      2
    else if (x == "SU") 
      3
    else 4
    vecFromList <- function(item, aList) {
      if (!is.list(aList[[1]])) 
        return(aList[[item]])
      else {
        tVec <- vector(length = 0)
        for (i in 1:length(aList)) {
          tVec <- append(tVec, (aList[[i]])[[item]])
        }
      }
      tVec
    }
    gamma <- vecFromList(1, parms)
    delta <- vecFromList(2, parms)
    xi <- vecFromList(3, parms)
    lambda <- vecFromList(4, parms)
    type <- vecFromList(5, parms)
    type <- sapply(type, tfun)
    N <- length(gamma)
    value <- .C("sJohnsonR", as.double(gamma), as.double(delta), 
                as.double(xi), as.double(lambda), as.integer(type), as.integer(N), 
                mn = double(N), med = double(N), mod = double(N), var = double(N), 
                third = double(N), fourth = double(N),PACKAGE="SuppDists")
    aList <- list(title = "Johnson Distribution", gamma = gamma, 
                  delta = delta, xi = xi, lambda = lambda, type = type)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
                 value$third, value$fourth, -1)
  }
