#' @title Cramer - von Mises statistics
#'
#' @description   Calculates the Cramer-von Mises test statistic \deqn{T(S_n)=\frac{1}{2q}\sum_{i=1}^{2q}\left(H^-_n(S_{n,i})-H^+_n(S_{n,i})\right)^2} where \eqn{H^-_n(\cdot)} and \eqn{H^+_n(\cdot)} are the empirical CDFs of the the sample of baseline covariates close to the cutoff from the left and right, respectively. See equation (12) in Canay and Kamat (2017).
#' @param S Numeric. The pooled sample of induced order statistics. The first column of S can be viewed as an independent sample of W conditional on Z being close to zero from the left. Similarly, the second column of S can be viewed as an independent sample of W conditional on Z being close to the cutoff from the right. See section 3 in Canay and Kamat (2017).
#' @return Returns the numeric value of the Cramer - von Mises test statistic.
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Canay, I and Kamat V. (2017) ``Approximate Permutation Tests and Induced Order Statistics in the Regression Discontinuity Design''
#' @keywords permutation test rdperm
#' @export


CvM.stat <- function(S){
  q <- nrow(S)
  tS <- apply(S,1,function(x) (H.cdf(S[,1],x[1])-H.cdf(S[,2],x[2]))^2)
  return((1/2*q)*sum(tS))
}
