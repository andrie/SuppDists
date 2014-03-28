#' A very long period pseudo-random generator.
#' 
#' A pseudo-random number generator with a period exceeding `10e9824` due to George Marsaglia
#' 
#' 
#' The following is an email from George Marsaglia on this generator:
#' 
#' `MCW1019` will provide random 32-bit integers at the rate of 300 million per second (on a 850MHz PC).
#' 
#' The period of MWC1029 exceeds 10e9824, making it billions and billions ... and billions times as long as the highly touted longest-period RNG, the Mersenne twister.  It is also several times as fast and takes a few lines rather than several pages of code.  (This is not to say that the Mersenne twister is not a good RNG; it is.  I just do not equate complexity of code with randomness.  It is the complexity of the underlying randomness that counts.)
#' 
#' As for randomness, it passes all tests in The Diehard Battery of Tests of Randomness \url{http://stat.fsu.edu/pub/diehard} as well as three new tough tests I have developed with the apparent property that a RNG that passes tuftsts.c will pass all the tests in Diehard.
#' 
#' MWC1019 has the property that every possible sequence of 1018 successive 32-bit integers will appear somewhere in the full period, for those concerned with the `equi-distribution` in dimensions `2,3,...1016,1017,1018`.
#' 
#' NOTE: This function does not work correctly on LP64 bit machines which use 64 bit longs.
#' 
#' The generator requires 1019 initial 32 bit random values to start. These are provided by using Marsaglia's MWC generator. By setting `new.start=FALSE`, the sequence may be sampled in blocks.  specifies an internal array dimension of 625 for seed length, and thus at the present time, it is not possible to implement this generator using the .Random.seed mechanism. In this implementation, `rMWC1019()` seems to be about twice as fast as the Mersenne-Twister implemented in \R.
#' 
#' @aliases rMWC1019
#' @param n number of values to generate. If n is a vector, length(n) values will be generated
#' @param new.start logical scalar. If TRUE the generator will be started afresh using the seed
#' @param seed scalar 32 bit integer used as the generator seed
#' @return generates a vector of random numbers.
#' @author Bob Wheeler \email{bwheelerg@@gmail.com}
#' @keywords distribution
#' @examples
#' 
#' 
#' rMWC1019(50,new.start=TRUE,seed=492166)
#' rMWC1019(50)
#' 
#' @export

rMWC1019 <- function (n, new.start = FALSE, seed = 556677) {
  n <- if (length(n) == 1) 
    n
  else length(n)
  .C("MWC1019R", val = double(n), as.integer(n), as.integer(new.start), 
     as.integer(seed),PACKAGE="SuppDists")$val
}
