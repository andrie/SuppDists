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
