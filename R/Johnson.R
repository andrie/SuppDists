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


pJohnson <- function (q, parms, lower.tail = TRUE, log.p = FALSE) 
{
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
