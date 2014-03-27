dghyper <- function (x, a, k, N, log = FALSE) {
  M <- max(length(x), length(a), length(k), length(N))
  x <- rep(x, length.out = M)
  a <- rep(a, length.out = M)
  k <- rep(k, length.out = M)
  N <- rep(N, length.out = M)
  value <- .C("dghyperR", as.integer(x), as.double(a), as.double(k), 
              as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}

pghyper <- function (q, a, k, N, lower.tail = TRUE, log.p = FALSE) {
  M <- max(length(q), length(a), length(k), length(N))
  q <- rep(q, length.out = M)
  a <- rep(a, length.out = M)
  k <- rep(k, length.out = M)
  N <- rep(N, length.out = M)
  if (lower.tail == TRUE) {
    value <- .C("pghyperR", as.integer(q), as.double(a), 
                as.double(k), as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }
  else {
    value <- .C("ughyperR", as.integer(q), as.double(a), 
                as.double(k), as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
  }
  if (log.p == TRUE) 
    value <- log(value)
  value
}

sghyper <- function (a, k, N) {
    M <- max(length(a), length(k), length(N))
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("sghyperR", as.double(a), as.double(k), as.double(N), 
                as.integer(M), mn = double(M), med = double(M), mod = double(M), 
                var = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    aList <- list(title = "Generalized Hypergeometric", a = a, 
                  k = k, N = N)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
                 value$third, value$fourth, -1)
  }


qghyper <- function (p, a, k, N, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(a), length(k), length(N))
    p <- rep(p, length.out = M)
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    value <- .C("qghyperR", as.double(p), as.double(a), as.double(k), 
                as.double(N), as.integer(M), val = double(M),PACKAGE="SuppDists")$val
    value
  }


rghyper <- function (n, a, k, N) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    K <- max(length(a), length(k), length(N))
    a <- rep(a, length.out = K)
    k <- rep(k, length.out = K)
    N <- rep(N, length.out = K)
    .C("rghyperR", as.double(a), as.double(k), as.double(N), 
       as.integer(n), as.integer(K), value = double(n),PACKAGE="SuppDists")$value
  }


tghyper <- function (a, k, N) {
    value <- .C("tghyperR", as.double(a), as.double(k), as.double(N), 
                strn =paste(rep(" ", 128), collapse=""),PACKAGE="SuppDists"  )
    value$strn
  }
