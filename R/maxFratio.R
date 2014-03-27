dmaxFratio <- function (x, df, k, log = FALSE) {
  if (log == TRUE) 
    p <- exp(p)
  N <- max(length(x), length(df), length(k))
  x <- rep(x, length.out = N)
  df <- rep(df, length.out = N)
  k <- rep(k, length.out = N)
  .C("dmaxFratioR", as.double(x), as.integer(df), as.integer(k), 
     as.integer(N), val = double(N),PACKAGE="SuppDists")$val
}


pmaxFratio <- function (q, df, k, lower.tail = TRUE, log.p = FALSE) {
    N <- max(length(q), length(df), length(k))
    q <- rep(q, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    if (lower.tail == TRUE) {
      value <- .C("pmaxFratioR", as.double(q), as.integer(df), 
                  as.integer(k), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    else {
      value <- .C("umaxFratioR", as.double(q), as.integer(df), 
                  as.integer(k), as.integer(N), val = double(N),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
      value <- log(value)
    value
  }


qmaxFratio <- function (p, df, k, lower.tail = TRUE, log.p = FALSE) {
    if (lower.tail == FALSE) 
      p <- 1 - p
    if (log.p == TRUE) 
      p <- exp(p)
    N <- max(length(p), length(df), length(k))
    p <- rep(p, length.out = N)
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    .C("qmaxFratioR", as.double(p), as.integer(df), as.integer(k), 
       as.integer(N), val = double(N),PACKAGE="SuppDists")$val
  }

rmaxFratio <- function (n, df, k) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- max(length(df), length(k))
    df <- rep(df, length.out = M)
    k <- rep(k, length.out = M)
    .C("rmaxFratioR", as.integer(df), as.integer(k), as.integer(n), 
       as.integer(M), value = double(n),PACKAGE="SuppDists")$value
  }

smaxFratio <- function (df, k) {
    N <- max(length(df), length(k))
    df <- rep(df, length.out = N)
    k <- rep(k, length.out = N)
    value <- .C("smaxFratioR", as.integer(df), as.integer(k), 
                as.integer(N), mn = double(N), med = double(N), mod = double(N), 
                var = double(N), third = double(N), fourth = double(N),PACKAGE="SuppDists")
    aList <- list(title = "Maximum F ratio", df = df, k = k)
    makeStatList(aList, value$mn, value$med, value$var, value$mod, 
                 value$third, value$fourth, 2)
  }
