#' Supplementary distributions and two random number generators.
#' 
#' 
#' The package includes ten distributions supplementing those built into \R:
#' 
#' \describe{ 
#' \item{\code{\link{dinvGauss}}}{The inverse Gaussian and Wald distributions}
#' \item{\code{\link{dKruskalWallis}}}{Kruskall-Wallis distribution}
#' \item{\code{\link{dKendall}}}{Kendall's Tau}
#' \item{\code{\link{dFriedman}}}{Friedman's chi squared}
#' \item{\code{\link{dSpearman}}}{Spearman's rho}
#' \item{\code{\link{dmaxFratio}}}{The maximum F-ratio distribution}
#' \item{\code{\link{dPearson}}}{Pearson product moment correlation coefficient}
#' \item{\code{\link{dJohnson}}}{Johnson distributions}
#' \item{\code{\link{dNormScore}}}{normal scores}
#' \item{\code{\link{dghyper}}}{Generalized hypergeometric distributions (Kemp and Kemp)}
#' }
#' 
#' In addition, two random number generators of George Marsaglia:
#' 
#' \describe{
#' \item{\code{\link{rMWC1019}}}{A very long period pseudo-random generator}
#' \item{\code{\link{rziggurat}}}{The Ziggurat normal and exponential generator}
#' }
#' 
#' @name SuppDists-package
#' @aliases SuppDists
#' @docType package
#' @useDynLib SuppDists
#' @author 
#' Original author: Bob Wheeler \email{bwheelerg@@gmail.com}; 
#' Maintainer: Andrie de Vries \email{apdevries@@gmail.com}
#' @keywords package
NULL
