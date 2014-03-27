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
