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
