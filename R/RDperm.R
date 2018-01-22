#' @title Regression Discontinuity Design Permutation Test
#'
#' @description A permutation test for continuity of covariates in Sharp Regression Discontinuity Design as described in Canay and Kamat (2017).
#'
#' @param W Character. Vector of covariates names. The procedure will test the null hypothesis of continuity of the distribution of each element in W at the cutoff.
#' @param z Character. Running variable name. This is the scalar random variable that defines, along with the cutoff, the treatment assignment rule in the sharp regression discontinuity design.
#' @param data Data.frame.
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. See remark 3.2 in Canay and Kamat (2017). The default is B=499.
#' @param q_type A fixed and small (relative to the sample size) natural number that will define the \eqn{q}{q} closest values of the order statistic of \eqn{Z}{Z} to the right and to the left of the cutoff. The default, 'rot', value is given by the feasible rule of thumb in footnote 4 of Canay and Kamat (2017), section 3.1. If 'arot', it calls for the Rule of Thumb described in equation (15) of Canay and Kamat (2017), section 3.1. The default option grows at a slower rate than the optional rule of thumb, but adds a larger constant.
#' @param cutoff Numeric. The scalar defining the threshold of the running variable.
#' @param test.statistic Character. A rank test statistic satisfying rank invariance. The default is a Cramer-von Mises test statistic.
#' @return The functions \code{summary} and \code{plot} are used to obtain and print a summary and plot of
#' the estimated regression discontinuity. The object of class \code{RDperm} is a list
#' containing the following components:
#'  \item{results}{Matrix. Test Statistic, P-values and Q}
#'  \item{test.statistic}{Test Statistic}
#'  \item{q_type}{Type of Q used in the calculations, can be either,  "Defined by User", the "Rule of Thumb"  or the "Alternative Rule of Thumb".}
#'  \item{n_perm}{number of permutations}
#'  \item{rv}{Character. Running variable name}
#'  \item{Z}{Vector. Running Variable}
#'  \item{cutoff}{cutoff}
#'  \item{data}{data set}
#'  \item{S}{Matrix. Pooled sample of induced order statistics}
#'  \item{S_perm}{List. Permutations of the induced order statistic.}
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Canay, I and Kamat V, (2017) Approximate Permutation Tests and Induced Order Statistics in the Regression Discontinuity Design. \url{http://faculty.wcas.northwestern.edu/~iac879/wp/RDDPermutations.pdf}
#' @keywords permutation test rdperm
#' @include H.cdf.R
#' @include CvM.stat.R
#' @import quantreg
#' @importFrom stats cor var runif
#' @examples
#' permtest<-RDperm(W=c("demshareprev"),z="difdemshare",data=lee2008)
#' summary(permtest)
#'\dontrun{
#' permtest<-RDperm(W=c("demshareprev","demwinprev"),z="difdemshare",data=lee2008)
#' summary(permtest)
#' }
#' @export


RDperm<-function(W,z,data,n.perm=499,q_type=10,cutoff=0,test.statistic="CvM"){
  W_z<-base::subset(data, select=c(W,z))
  colnames(W_z)[colnames(W_z)==z]<-"z"
  N<-dim(data)[1]

  W_left <- W_z[W_z$z<cutoff,]
  n_left <- length(W_left$z)
  W_right <- W_z[W_z$z>=cutoff,]
  n_right <- length(W_right$z)

  if(N!=n_left+n_right) stop( paste( "Something is wrong with the number of observations", sep="" ) )

  # Induced order of W obs
  W_left <- W_left[order(W_left$z),]
  W_right <- W_right[order(W_right$z),]

  
  
  if(!(q_type%in%c("rot","arot")) & length(W)>1 ) {
    results<-matrix(NA, nrow=length(W)+1, ncol=3)
  } else results<-matrix(NA, nrow=length(W), ncol=3)


  # Selecting Q,
  if(q_type%in%c("rot","arot")){
    w<-as.list(W)
    if(q_type=="rot") rot<-lapply(w,qrot,W_z)
    if(q_type=="arot") rot<-lapply(w,aqrot,W_z)
    w<-mapply(c, w, rot, SIMPLIFY=FALSE)

    test<-lapply(w,function(x) {
      f<-RDperm.base(x[1],W_left, n_left, W_right, q=as.numeric(x[2]), z, n.perm, test.statistic)
      ret<-list()
      ret$test_statistic.obs<-f$test_statistic.obs
      ret$pvalues<-f$pvalues
      ret$q<-as.numeric(f$q)
      return(ret)
    })
    test<-do.call(rbind,test)
    test<-t(apply(test,1,function(x) do.call(rbind,x)))

    results[,1]<-test[,1]
    results[,2]<-test[,2]
    results[,3]<-test[,3]

  }

  if(!(q_type%in%c("rot","arot")) & length(W)>1 )  results[,3]<-rep(q_type,length(W)+1)
  if(!(q_type%in%c("rot","arot")) & length(W)==1 ) results[,3]<-rep(q_type,length(W))
  q<-min(results[,3])

  permtest<-RDperm.base(W,W_left, n_left, W_right, q=q,z, n.perm, test.statistic)

  if(!(q_type%in%c("rot","arot")) ){
  results[,1]<-permtest$test_statistic.obs
  results[,2]<-permtest$pvalues
  }

  if((q_type%in%c("rot","arot")) & length(W)>1 ){
    permtest<-RDperm.base(W,W_left, n_left, W_right, q=q,z, n.perm, test.statistic)
    results_updated<-matrix(NA, nrow=length(W)+1, ncol=3)
    results_updated[,1]<-c(results[,1],permtest$test_statistic.obs[length(W)+1])
    results_updated[,2]<-c(results[,2],permtest$pvalues[length(W)+1])
    results_updated[,3]<-c(results[,3],q)
    results<-results_updated
  }

  for(i in 1:3) results[,i]<- as.numeric(results[,i])
  colnames(results)<-c("T(Sn)","Pr(>|z|)", "q")
  if(length(W)>1){rownames(results)<-c(W,"Joint.Test")
  }else rownames(results)<-W


  object_perm<-list()


  object_perm$results<-results
  object_perm$test.statistic<-test.statistic

  object_perm$Z<- W_z[,"z"] #running variable
  object_perm$rv<- z #name of running variable
  object_perm$cutoff<- cutoff #cutoff


  if(q_type=="rot"){
    object_perm$q_type<-"Rule of Thumb"
  } 
  else if(q_type=="arot"){
    object_perm$q_type<-"Alternative Rule of Thumb"
  }
    else{object_perm$q_type<- "Defined by User"}

  object_perm$n_perm<- n.perm #number of permutations
  object_perm$data<-data


  class(object_perm)<-"RDperm"

  return(object_perm)
}



RDperm.base<-function(W,W_left, n_left, W_right, z, q, n.perm, test.statistic){

  #Step 1 & 2. Compute the order statistics of Z and the associated values of W
  # q closest W from the left/right of threshold
  # Equation (10)
  W_left_q<-base::subset(W_left[(n_left-q+1):n_left,], select=c(W))
  Z_left<-base::subset(W_left[(n_left-q+1):n_left,], select=c(z))
  W_right_q<-base::subset(W_right[1:q,], select=c(W))
  Z_right <-base::subset(W_right[1:q,], select=c(z))

  Sn<-rbind(W_left_q,W_right_q)



  if(test.statistic=="CvM"){
    #Step 3. Compute the test statistic
    test_statistic.obs<-apply(Sn,2,CvM.stat)
    if(length(W)>1){
    n.test_statistic.obs<-names(test_statistic.obs)
    K<-length(W)
    c<-C.unitsphere(K)
    cS<-as.matrix(Sn)%*%c
    TSn.joint<-max(apply(cS,2,calc_stat.CvM))
    test_statistic.obs<-c(test_statistic.obs,TSn.joint)
    names(test_statistic.obs)<-c(n.test_statistic.obs,"joint")
    }

    #Step 4. Generate random permutations
    sample.indexes = lapply(1:n.perm, function(x) sample(1:(2*q)))
    S_perm_list<-lapply(sample.indexes,function(x,db) {db[x,]},Sn)

    calc_stat_res<-lapply(S_perm_list,calc_stat.CvM)


    #Step 6. Compute the p-value of the test
    test_statistic<-"CvM"
    ind.rule<-lapply(calc_stat_res,function(x,CvM) {ifelse(x>=CvM,1,0)},test_statistic.obs)
    ind.rule<-do.call(cbind,ind.rule)
    ind.rule<-rowMeans(ind.rule)
  } else{"Need to generate Kologorov Statistic"}

  object_perm<-list()
  object_perm$test_statistic.obs<-test_statistic.obs
  object_perm$pvalues<-ind.rule
  object_perm$q<- q #q
  #object_perm$S<- S_perm
  object_perm$S<- Sn
  object_perm$S_perm<- S_perm_list


  class(object_perm)<-"RDperm"

  return(object_perm)
}





calc_stat.CvM<-function(x){
  if(is.vector(x)==T){
    stat<-CvM.stat(x)
  }
  else {
    stat<-apply(x,2,CvM.stat)
    n.stat<-names(stat)
    K<-dim(x)[2]
    c<-C.unitsphere(K)
    cS<-as.matrix(x)%*%c
    TSn.joint<-max(apply(cS,2,calc_stat.CvM))
    stat<-c(stat,TSn.joint)
    names(stat)<-c(n.stat,"joint")
    }
  return(stat)
}



aqrot<-function(w,W_z){
  w<-W_z[,w]
  z<-W_z[,"z"]
  N<-length(w)
  t <- seq(from=min(z),to=max(z),length.out = 2*N)
  f <- quantreg::akj(z,t)$dens
  t0 <- which(abs(t-0)==min(abs(t-0)))
  f0 <- f[t0]
  q<-ceiling(f0*var(z)*sqrt((1-cor(z,w)^2))*(N^{0.9}/log(N)))
  if(q<10){
    q<-10
  }else if(q>N^(0.9)/log(N)){
    q<-ceiling(N^(0.9)/log(N))
  } else {
    q<-q
  }
}



qrot<-function(w,W_z){
  w<-W_z[,w]
  z<-W_z[,"z"]
  N<-length(w)
  t <- seq(from=min(z),to=max(z),length.out = 2*N)
  f <- quantreg::akj(z,t)$dens
  t0 <- which(abs(t-0)==min(abs(t-0)))
  f0 <- f[t0]
  q<-ceiling(f0*var(z)*sqrt(10*(1-cor(z,w)^2))*(N^{0.75}/log(N)))
  if(q<10){
    q<-10
  }else if(q>N^(0.9)/log(N)){
    q<-ceiling(N^(0.9)/log(N))
  } else {
    q<-q
  }
}


C.unitsphere <- function(K){
  # Store Matrix
  C <- matrix(NA,nrow = K,ncol=100-K)
  # Fill the matrix with random numbers in [-1,1]
  C <- apply(C,2,function(x) stats::runif(x,-1,1))
  # Normalize it so each column c of C is such that ||c||=1
  C <- apply(C,2,function(x) x/sum(x))
  # Return C and the K canonical elemnts (vectors with zeros in
  # all coordinates except for one)
  return(cbind(C,diag(K)))
}
