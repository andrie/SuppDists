#' Internal SuppDists function.
#' 
#' This is not to be called by the user.
#' 
#' @aliases makeStatList
#' @keywords distribution internal
makeStatList <- function (head, mn, med, var, mod, third, fourth, dig) {
    sd <- sqrt(var)
    skew <- sign(third) * abs(third)/sd^3
    kurt <- -3 + fourth/var^2
    pskew <- (mn - mod)/sd
    if (dig > 0) {
      mn <- round(mn, digits = dig)
      med <- round(med, digits = dig)
      mod <- round(mod, digits = dig)
      var <- round(var, digits = dig)
      sd <- round(sd, digits = dig)
      third <- round(third, digits = dig)
      fourth <- round(fourth, digits = dig)
      pskew <- round(pskew, digits = dig)
      skew <- round(skew, digits = dig)
      kurt <- round(kurt, digits = dig)
    }
    theList <- list(Mean = mn, Median = med, Mode = mod, Variance = var, 
                    SD = sd, ThirdCentralMoment = third, FourthCentralMoment = fourth, 
                    PearsonsSkewness...mean.minus.mode.div.SD = pskew, Skewness...sqrtB1 = skew, 
                    Kurtosis...B2.minus.3 = kurt)
    c(head, theList)
  }

moments <- function (x) {
    N <- length(x)
    v <- ((N - 1)/N) * var(x)
    sigma <- sqrt(v)
    m3 <- (sum((x - mean(x))^3))/N
    skew <- m3/sigma^3
    m4 <- (sum((x - mean(x))^4))/N
    kurt <- (m4/v^2) - 3
    c(mean = mean(x), sigma = sigma, skew = skew, kurt = kurt)
  }

normOrder <- function (N) {
    N <- if (length(N) > 1) 
      length(N)
    else N
    M <- N%/%2
    value <- .C("normOrdR", val = double(M), as.integer(N), as.integer(M),PACKAGE="SuppDists")$val
    if (0 == N%%2) 
      c(-value, rev(value))
    else c(-value, 0, rev(value))
  }




