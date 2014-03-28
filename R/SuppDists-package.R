#' Supplementary distributions and two random number generators.
#' 
#' 
#' The package includes ten distributions supplementing those built into \R:
#' 
#' \describe{ 
#' \item{\code{invGauss}}{The inverse Gaussian and Wald distributions}
#' \item{\code{KruskalWallis}}{Kruskall-Wallis distribution}
#' \item{\code{Kendall}}{Kendall's Tau}
#' \item{\code{Friedman}}{Friedman's chi squared}
#' \item{\code{Spearman}}{Spearman's rho}
#' \item{\code{maxFratio}}{The maximum F-ratio distribution}
#' \item{\code{Pearson}}{Pearson product moment correlation coefficient}
#' \item{\code{Johnson}}{Johnson distributions}
#' \item{\code{NormalScore}}{normal scores}
#' \item{\code{ghyper}}{Generalized hypergeometric distributions (Kemp and Kemp)}
#' }
#' 
#' In addition, two random number generators of George Marsaglia:
#' 
#' \describe{
#' \item{\code{MWC1019}}{A very long period pseudo-random generator}
#' \item{\code{ziggurat}}{The Ziggurat normal and exponential generator}
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
