#' @title Regression Discontinuity Design Permutation test
#'
#' @description  Calculates the empirical CDF of the sample of \eqn{W}{W} conditional on \eqn{Z}{Z} being close to the cutoff from either the left or right. Given the induced order for the baseline covariates \deqn{W^{-}_{[q]}, W^{-}_{[q-1]},\dots\le W^{-}_{[1]}} or \deqn{W^{+}_{[1]}, W^{+}_{[2]},\dots, W^{+}_{[q]}}, this function will calculate either \deqn{H^-_n(t)=\frac{1}{q}\sum_{i=1}^q I\{W^{-}_{[i]}\le t\}}  or  \deqn{H^+_n(t)=\frac{1}{q}\sum_{i=1}^q I\{W^{+}_{[i]}\le t\}} depending on the argument of the function. See section 3 in Canay & Kamat (2017).
#' @param W Numeric. The sample of induced order statistics. The input can be either \eqn{\{W^{-}_{[q]}, W^{-}_{[q-1]},\dots, W^{-}_{[1]}\}} or \eqn{\{W^{+}_{[1]}, W^{+}_{[2]},\dots, W^{+}_{[q]}\}}.
#' @param t Numeric. The scalar needed for the calculation of the CDF.
#' @return Numeric. For a sample \eqn{W=(w_1,\dots,w_n)}, returns the fraction of observations less or equal to \eqn{t}{t}.
#'
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Canay, I and Kamat V, (2017) Approximate Permutation Tests and Induced Order Statistics in the Regression Discontinuity Design. \url{http://faculty.wcas.northwestern.edu/~iac879/wp/RDDPermutations.pdf}
#' @keywords permutation test rdperm
#' @export

H.cdf <- function(W,t){
  q <- length(W)
  H <- (1/q)*sum((W<=t))
  return(H)
}
