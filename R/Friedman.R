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
