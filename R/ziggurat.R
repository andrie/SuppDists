#' The Ziggurat normal and exponential generator.
#' 
#' Generates normal and exponential random pseudo-random numbers by the method of Marsaglia and Tsang.
#' 
#' 
#' @aliases rziggurat
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param normal logical scalar; if TRUE normal values are produced, otherwise exponential values
#' @param new.start logical scalar. If TRUE the generator will be started afresh using the seed
#' @param seed scalar 32 bit starting integer
#' @return Generates a vector of real pseudo random numbers.
#' @note This function does not work properly on LP64 machines which use 64 bit longs.
#' 
#' This implementation running in is approximately three times as fast as `rnorm()`.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @references
#' Marsaglia, George, and Tsang, Wai Wan. 2000. The Ziggurat method for generating random variables. 
#' \emph{Journal of Statistical software.} \bold{5-8.} \url{http://www.jstatsoft.org/v05/i08/}
#' @keywords distribution
#' @examples
#' 
#' rziggurat(50, new.start=TRUE)
#' rziggurat(50)
#' rziggurat(50, new.start=TRUE)
#' 
#' @export

rziggurat <- function (n, normal = TRUE, new.start = FALSE, seed = 556677) {
  n <- if (length(n) > 1) 
    length(n)
  else n
  .C("ziggR", val = double(n), as.integer(n), as.integer(normal), 
     as.integer(new.start), as.integer(seed),PACKAGE="SuppDists")$val
}
