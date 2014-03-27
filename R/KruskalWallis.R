dKruskalWallis <- function (x, c, N, U, log = FALSE) {
  M <- max(length(x), length(c), length(N), length(U))
  x <- rep(x, length.out = M)
  c <- rep(c, length.out = M)
  n <- rep(N, length.out = M)
  U <- rep(U, length.out = M)
  Ns <- rep(FALSE, length.out = M)
  value <- .C("dKruskalWallisR", as.double(x), as.integer(c), 
              as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
              val = double(M),PACKAGE="SuppDists")$val
  if (log == TRUE) 
    value <- log(value)
  value
}


pKruskalWallis <- function (q, c, N, U, lower.tail = TRUE, log.p = FALSE) {
    M <- max(length(q), length(c), length(N), length(U))
    q <- rep(q, length.out = M)
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    if (lower.tail == TRUE) {
      value <- .C("pKruskalWallisR", as.double(q), as.integer(c), 
                  as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
                  val = double(M),PACKAGE="SuppDists")$val
    }
    else {
      value <- .C("uKruskalWallisR", as.double(q), as.integer(c), 
                  as.integer(n), as.double(U), as.integer(Ns), as.integer(M), 
                  val = double(M),PACKAGE="SuppDists")$val
    }
    if (log.p == TRUE) 
      value <- log(value)
    value
  }


qKruskalWallis <- function (p, c, N, U, lower.tail = TRUE, log.p = FALSE) {
    if (log.p == TRUE) 
      p <- exp(p)
    if (lower.tail == FALSE) 
      p <- 1 - p
    M <- max(length(p), length(c), length(N), length(U))
    p <- rep(p, length.out = M)
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    .C("qKruskalWallisR", as.double(p), as.integer(c), as.integer(N), 
       as.double(U), as.integer(Ns), as.integer(M), val = double(M),
       PACKAGE="SuppDists")$val
  }


rKruskalWallis <- function (n, c, N, U) {
    n <- if (length(n) > 1) 
      length(n)
    else n
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    N <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    .C("rKruskalWallisR", randArray = double(n), as.integer(n), 
       as.integer(M), as.integer(c), as.integer(N), as.double(U), 
       as.integer(Ns),PACKAGE="SuppDists" )$randArray
  }


sKruskalWallis <- function (c, N, U) {
    M <- max(length(c), length(N), length(U))
    c <- rep(c, length.out = M)
    n <- rep(N, length.out = M)
    U <- rep(U, length.out = M)
    Ns <- rep(FALSE, length.out = M)
    value <- .C("sKruskalWallisR", as.integer(c), as.integer(n), 
                as.double(U), as.integer(Ns), as.integer(M), var = double(M), 
                mod = double(M), third = double(M), fourth = double(M),PACKAGE="SuppDists")
    mn <- (c - 1)
    aList <- list(title = "Kruskal Wallis", c = c, N = n, U = U)
    median <- qKruskalWallis(0.5, c, n, U, Ns)
    makeStatList(aList, mn, median, value$var, value$mod, value$third, 
                 value$fourth, -1)
  }
