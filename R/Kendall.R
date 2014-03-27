dKendall <- function (x, N, log = FALSE) {
  M <- max(length(x), length(N))
  x <- rep(x, length.out = M)
  N <- rep(N, length.out = M)
  value <- .C("dKendallR", as.integer(N), as.double(x), as.integer(M), 
              val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}


pKendall <- function (q, N, lower.tail = TRUE, log.p = FALSE) {
  M <- max(length(q), length(N))
  q <- rep(q, length.out = M)
  N <- rep(N, length.out = M)
  if (lower.tail == TRUE) {
    value <- .C("pKendallR", as.integer(N), as.double(q), 
                as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }
  else {
    value <- .C("uKendallR", as.integer(N), as.double(q), 
                as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }
  if (log.p == TRUE) 
    value <- log(value)
  value
}


qKendall <- function (p, N, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(N))
    p <- rep(p, length.out = M)
    N <- rep(N, length.out = M)
    .C("qKendallR", as.integer(N), as.double(p), as.integer(M), 
       val = double(M),PACKAGE="SuppDists")$val
  }

rKendall <- function (n, N) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- length(N)
    .C("rKendallR", as.integer(N), as.integer(n), as.integer(M), 
       val = double(n),PACKAGE="SuppDists")$val
  }

sKendall <- function (N) {
    M <- length(N)
    mn <- rep(0, length.out = M)
    med <- rep(0, length.out = M)
    mod <- rep(0, length.out = M)
    third <- rep(0, length.out = M)
    var <- (4 * N + 10)/(9 * N * (N - 1))
    fourth <- .C("fourthKendallR", as.integer(N), as.integer(M), 
                 val = double(M),PACKAGE="SuppDists")$val
    aList <- list(title = "Kendall's Tau", N = N)
    makeStatList(aList, mn, med, var, mod, third, fourth, -1)
  }
