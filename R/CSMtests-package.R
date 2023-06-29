#' CSMtests
#'
#' Tools for central statistical monitoring.
#'
#' At the moment nothing but correlation-based tests.
#'
#' @seealso \code{\link{testCorFisher}}, \code{\link{testCorMargin}}m \code{\link{calcPars}}
#'
#'@examples
#'Np = round(runif(100)*50+10);
#'x = simulateCentersCors(.5, .1, Np)
#'ctr = rep(1:100, Np);
#'
#'pp = calcPars(x, ctr)
#'
#'p1 = testCorFisher(pp$cc, pp$N)
#'p2 = testCorMargin(pp$cc, pp$N, pp$sd)
#' @docType package
#' @name CSMtests
NULL

#' @import MASS
NULL

#' Distribution of center sizes in real studies.
#'
#' Distribution of center sizes (N patients) for 3 real studies of different sizes.
#' 
#'
#' @format ## `centerSizes`
#' A list of length 3 (small, medium, large). 
#' Each item of the list is the N patients of a center from that study
"centerSizes"