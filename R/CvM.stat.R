#' @title Cramer - von Mises statistics
#'
#' @description   Calculates the Cramer-von Mises test statistic \deqn{T(S_n)=\frac{1}{2q}\sum_{i=1}^{2q}\left(H^-_n(S_{n,i})-H^+_n(S_{n,i})\right)^2} where \eqn{H^-_n(\cdot)} and \eqn{H^+_n(\cdot)} are the empirical CDFs of the the sample of baseline covariates close to the cutoff from the left and right, respectively. See equation (12) in Canay and Kamat (2017).
#' @param Sn Numeric. The pooled sample of induced order statistics. The first column of S can be viewed as an independent sample of W conditional on Z being close to zero from the left. Similarly, the second column of S can be viewed as an independent sample of W conditional on Z being close to the cutoff from the right. See section 3 in Canay and Kamat (2017).
#' @return Returns the numeric value of the Cramer - von Mises test statistic.
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Canay, I and Kamat V, (2017) Approximate Permutation Tests and Induced Order Statistics in the Regression Discontinuity Design. \url{http://faculty.wcas.northwestern.edu/~iac879/wp/RDDPermutations.pdf}
#' @keywords permutation test rdperm
#' @export



CvM.stat <- function(Sn){
  q <- length(Sn)/2

  Ind_left<-outer(Sn[1:q],Sn,"<=")
  H_left<-apply(Ind_left,2,sum)/q

  Ind_right<-outer(Sn[(q+1):(2*q)],Sn,"<=")
  H_right<-apply(Ind_right,2,sum)/q

  return(sum((H_left-H_right)^2)/(2*q))
}

