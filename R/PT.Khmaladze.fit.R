
#' @title  Non-Parametric Hypothesis Testing with a Nuisance Parameter: A Permutation Test
#'
#' @description A permutation test of the two-sample goodness-of-fit hypothesis in the presence of an estimated niusance parameter. The permutation test considered here is based on the Khmaladze transformation of the empirical process (Khmaladze (1981)), and adapted by Chung and Olivares-Gonzalez (2018).
#' 
#' @param y1 Numeric. A vector containing the response variable of the treatment group.
#' @param y0 Numeric. A vector containing the response variable of the control group.  
#' @param alpha Numeric. Nominal level for the test. The default is 0.05.
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. The default is n.perm=999.

#' @return An object of class "PT.Khmaladze.fit" is a list containing at least the following components:
#' 
#'  \item{n_populations}{Number of grups.}
#'  \item{N}{Sample Size.}
#'  \item{T.obs}{Observed test statistic.}
#'  \item{shift}{The estimated nuisance parameter (average treatment effect).}
#'  \item{cv}{Critical Value. This value is used in the general construction of a randomization test.}
#'  \item{pvalue}{P-value.}
#'  \item{T.perm}{Vector. Test statistic recalculated for all permutations used in the stochastic approximation.}
#'  \item{n_perm}{Number of permutations.}
#'  \item{sample_sizes}{Groups size.}
#'  
#' @author Maurcio Olivares-Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references  
#' Khmaladze, E. (1981). Martingale Approach in the Theory of Goodness-of-fit Tests. Theory of Probability and its Application, 26: 240–257.
#' Chung, Eunyi and Mauricio Olivares (2018). Non-Parametric Hypothesis Testing with a Nuisance Parameter: A Permutation Test Approach. Working Paper.
#' 
#' @keywords permutation test goodness-of-fit Khmaladze Transformation
#' @import quantreg
#' @importFrom stats quantile cor var runif 
#' @examples
#'\dontrun{
#' Y0 <- rnorm(100, 1, 1)
#' # Treatment Group with constant shift equals to 1
#' Y1 <- Y0 + 1
#' Tx = sample(100) <= 0.5*(100)
#' # Observed Outcome 
#' Y = ifelse( Tx, Y1, Y0 )
#' dta <- data.frame(Y = Y, Z = as.numeric(Tx))
#' pt.GoF<-PT.Khmaladze.fit(dta$Y[dta$Z==1],data$Y[dta$Z==0],n.perm = 49)
#' summary(pt.GoF)
#' }
#' @export



PT.Khmaladze.fit <- function(y1, y0, alpha=0.05, n.perm=999){

  if (anyNA(c(y1,y0))) stop("NAs in first or second argument")
  if (!is.numeric(c(y1,y0)))  stop("Arguments must be numeric")
  
  # ---------------------------- #
  #   Calculates the statistic   #
  # ---------------------------- #

  # Sample Size
  N<-length(c(y1,y0)) 
  m <- length(y1)
  lengths <- list(m,N-m)
  # Calculate the shift (average treatment effect)
  shift <- mean(y1)-mean(y0) 
  # Recentering
  Y1.star <- y1- shift
  Y0.star <- y0
  # Grid needed for the density estimation and numeric integration
  p <- 3*N
  eps <- .Machine$double.eps
  t <- seq(eps, 1, length.out=p)
  # Permutations of data
  Sn <- c(Y1.star,Y0.star)
  sample.indexes = lapply(1:n.perm, function(x) sample(1:N))
  S_perm_list<-lapply(sample.indexes,function(x,db) {db[x]},Sn)
  S_perm<-do.call(cbind,S_perm_list)
  Sn <- cbind(Sn,S_perm)
    
  # Calculation of the test statistic for all permutations. 
  stat <- apply(Sn,2, function(x) Khm.trans(x,m,N-m,n.perm,t,p))
  
  # Observed Test statistic
  T.obs<-stat[1] 
  #Test statistic for the permutated samples
  T.perm<-stat[-1] 
  
  # Critical Value
  cv <- randomization.test(stat,alpha)
  #Indicator rule
  ind.rule<- mean(ifelse(T.perm>=T.obs,1,0))
  
  
  object_perm<-list() #Generates an empty list to collect all the required info for summary
  object_perm$N<-N
  object_perm$T.obs<-T.obs
  object_perm$shift <- shift
  object_perm$cv <- cv
  object_perm$pvalue<-ind.rule
  object_perm$T.perm<- T.perm
  object_perm$n_perm<- n.perm
  object_perm$sample_sizes<-c(lengths)
  
  
  class(object_perm)<-"PT.Khmaladze.fit"
  
  return(object_perm)
  
}

"Khm.trans" <- function(Z,m,n,n.perm,t,p){
  # Define treatent and control groups (these are recentered values)
  Y1 <- Z[1:m]
  Y0 <- Z[(m+1):(m+n)]
  # Extended score, nonparametric. Calls for quantreg function akj
  score <- akj(Y0,t)
  gdot2 <- -score$psi*score$dens
  gdot0 <- cbind(rep(1, (p-1)), gdot2[1:(length(gdot2)-1)])
  # Transformations needed for the calculation of the uniform empirical process
  Y1.s <- ecdf(Y1)
  F1 <- Y1.s(Y1)
  Y0.s <- ecdf(Y0)
  F0 <- Y0.s(Y0)
  g1 <- ecdf(quantile(Y0,t))(t)
  Vn.hat1 <- ecdf(F1)
  Vn.hat0 <- ecdf(F0)
  # Calculation of the compensator 
  comp <- K.hat(t,p,Y1,Y0,gdot0)
  
  # The Khmaladze transformation of the uniform empirical process
  # Intuitively, this is like calculating the residuals in a regression
  # of the uniform empirical process on the (extended) score function.
  return(sqrt(n*m/(m+n))*max(abs((Vn.hat1(g1) -Vn.hat0(t)-comp ))))
}

"K.hat" <- function(t, p, Y1, Y0, gdot) {
  # Calculates the compensator using recursive least squares
  # The compensator is like the fitted values in a regression of the 
  # uniform empirical process on the (extended) score function.
  
  # The "covariates"
  dt <- diff(t)
  X <- gdot#/sqrt(dt)
  X <- X[(p-1):1, ]
  
  # The "dependent variable"
  Y1.s <- ecdf(Y1)
  F1 <- Y1.s(Y1)
  Y0.s <- ecdf(Y0)
  F0 <- Y0.s(Y0)
  f1 <- ecdf(quantile(Y0,t))(t)
  Vn.hat1 <- ecdf(F1)
  Vn.hat0 <- ecdf(F0)
  y <- rev((diff(Vn.hat1(f1)-Vn.hat0(t))))#*sqrt(dt)
  
  # Calculate the regression coefficients using recursive least squares
  bhat <- lm.fit.recursive(X,y, int=FALSE)
  bhat <- bhat[ , (p-1):1]
  # Fitted values
  dk <- diag(gdot%*%bhat)
  comp <- c(0,cumsum(dk))
  return(comp)
}



"randomization.test" <- function(x,alpha){
  # Number of permutations
  n.perm <- length(x)
  # Observed test statistic
  Tn.obs <- x[1]
  
  # Corret size of the test
  # See Lehmann & Romano (2005), p. 633
  y <- sort(x)
  k <- n.perm - floor(n.perm*alpha)
  cv = y[k]
  M.plus <- sum(y>cv)
  M.zero <- sum(y==cv)
  a <-  (n.perm*alpha - M.plus)/M.zero
  
  # Rejection
  if(Tn.obs > cv){
    out <-  1
  } 
  else if(Tn.obs == cv){
    out <-  ifelse(runif(1) <= a, 1, 0)
  } 
  else {
    out <-  0 
  }
  return(c(cv,out))
}











