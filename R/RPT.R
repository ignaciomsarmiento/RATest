#' @title Robust Permutation Test
#'
#' @description This function considers the general problem of inference from the permutation distribution when comparing parameters from k populations. The test statistics will be based on the difference of estimators that are asymptotically linear. For illustrative purposes we will consider here the 2 sample case, but the function works for k-samples. 
#'
#' Difference of means: Here, the null hypothesis is of the form \eqn{H_0: \mu(P)-\mu(Q)=0}, and the corresponding test statistic is given by 
#'  \deqn{T_{m,n}=\frac{N^{1/2}(\bar{X}_m-\bar{Y}_n)}{\sqrt{\frac{N}{m}\sigma^2_m(X_1,\dots,X_m)+ \frac{N}{n}\sigma^2_n(Y_1,\dots,Y_n)}}}  
#'  where \eqn{\bar{X}_m} and \eqn{\bar{Y}_n} are the sample means from population \eqn{P} and population \eqn{Q}, respectively, and \eqn{\sigma^2_m(X_1,\dots,X_m)} is a consistent estimator of \eqn{\sigma^2(P)$ when $X_1,\dots,X_m} are i.i.d. from \eqn{P}. Assume consitency also under \eqn{Q}.
#'
#' Difference of medians: Let \eqn{F} and \eqn{G} be the CDFs corresponding to \eqn{P} and \eqn{Q}, and denote \eqn{\theta(F)} the median of \eqn{F} i.e. \eqn{\theta(F)=\inf\{x:F(x)\ge1/2\}}. Assume that \eqn{F} is continuously differentiable at \eqn{\theta(P)} with derivative \eqn{F'} (and the same with \eqn{F} replaced by \eqn{G}). Here, the null hypothesis is of the form  \eqn{H_0: \theta(P)-\theta(Q)=0}, and the corresponding test statistic is given by
#' \deqn{T_{m,n}=\frac{N^{1/2}\left(\theta(\hat{P}_m)-\theta(\hat{Q})\right)}{\hat{\upsilon}_{m,n}}}
#' where \eqn{\hat{\upsilon}_{m,n}} is a consistent estimator of \eqn{\upsilon(P,Q)}:
#' \deqn{\upsilon(P,Q)=\frac{1}{\lambda}\frac{1}{4(F'(\theta))^2}+\frac{1}{1-\lambda}\frac{1}{4(G'(\theta))^2}}
#' Choices of \eqn{\hat{\upsilon}_{m,n}} may include the kernel estimator of Devroye and Wagner (1980), the bootstrap estimator of Efron (1992), or the smoothed bootstrap  Hall et al. (1989) to list a few. For further details, see Chung and Romano (2013). Current implementation uses the bootstrap estimator of Efron (1992)
#'                                                                                        
#'Difference of variances: Here, the null hypothesis is of the form \eqn{H_0: \sigma^2(P)-\sigma^2(Q)=0}, and the corresponding test statistic is given by 
#'\deqn{T_{m,n}=\frac{N^{1/2}(\hat{\sigma}_m^2(X_1,\dots,X_,)-\hat{\sigma}_n^2(Y_1,\dots,Y_n))}{\sqrt{\frac{N}{m}(\hat{\mu}_{4,x}-\frac{(m-3)}{(m-1)}(\hat{\sigma}_m^2)^2)+\frac{N}{n}(\hat{\mu}_{4,y}-\frac{(n-3)}{(n-1)}(\hat{\sigma}_y^2)^2)}}} 
#'where \eqn{\hat{\mu}_{4,m}} the sample analog of \eqn{E(X-\mu)^4} based on an iid sample \eqn{X_1,\dots,X_m} from \eqn{P}. Similarly for \eqn{\hat{\mu}_{4,n}}.
#'
#'
#' @param formula a formula object, in which the response variable is on the left of a ~ operator, and the groups on the right.
#' @param data a data.frame containing the named variables needed for the formula. If this argument is missing, then the variables in the formula should be on the search list. 
#' @param test testing problem. It admits  "means" if the objective is to test for difference of Means, "medians" for difference of Medians, and "variances" for difference of Variances. In the case the user is interested in testing for difference of medians, the Efron (1992) bootstrap estimator is used to estimate the variances (For further details, see Chung and Romano (2013))
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. The default is n.perm=499.	
#' @param na.action a function to filter missing data. This is applied to the model.frame . The default is na.omit, which deletes observations that contain one or more missing values.

#' @return An object of class "RPT" is a list containing at least the following components:


#'  \item{description}{Type of test, can be Difference of Means, Medians, or Variances.}
#'  \item{n_populations}{Number of grups.}
#'  \item{N}{Sample Size.}
#'  \item{T.obs}{Observed test statistic.}
#'  \item{pvalue}{P-value.}
#'  \item{T.perm}{Vector. Test statistics calculated from the permutations of the data.}
#'  \item{n_perm}{Number of permutations.}
#'  \item{parameters}{Estimated parameters.}
#'  \item{sample_sizes}{Sample size of groups.}

#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references 
#' Chung, E. and Romano, J. P. (2013). Exact and asymptotically robust permutation tests. The Annals of Statistics, 41(2):484–507.
#' 
#' Chung, E. and Romano, J. P. (2016). Asymptotically valid and exact permutation tests based on two-sample u-statistics. Journal of Statistical Planning and Inference, 168:97–105.
#'
#' Devroye, L. P. and Wagner, T. J. (1980). The strong uniform consistency of kernel density estimates. In Multivariate Analysis V: Proceedings of the fifth International Symposium on Multivariate Analysis, volume 5, pages 59–77.
#' 
#' Efron, B. (1992). Bootstrap methods: another look at the jackknife. In Breakthroughs in statistics, pages 569–593. Springer.
#' Hall, P., DiCiccio, T. J., and Romano, J. P. (1989). On smoothing and the bootstrap. The Annals of Statistics, pages 692–704.

#' @keywords robust permutation test rpt
#' @import quantreg
#' @importFrom stats median pbinom
#' @examples
#'\dontrun{
#' male<-rnorm(50,1,1)
#' female<-rnorm(50,1,2)
#' dta<-data.frame(group=c(rep(1,50),rep(2,50)),outcome=c(male,female))
#' rpt.var<-RPT(dta$outcome~dta$group,test="variances")
#' summary(rpt.var)
#' 
#' }
#' @export



RPT<-function(formula,data,test="means",n.perm=499, na.action){
  if(!(test%in%c("means","medians","variances"))){
    print("Must specify test to perform, either 'means', 'medians', or 'variances'")
    stop()
    }
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data",  "na.action"), names(mf), 0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
 
  mf[,2] <- as.factor(mf[,2])
    
  N<-dim(mf)[1]
  
  k<-length(levels(mf[,2]))
  
  names<-as.character(unique(mf[,2]))
  
  #Permutations
  Sn<-mf[,1]
  groups<-mf[,2]
  mf<-mf[,c(2,1)]
  
  # Generate random permutations
  sample.indexes = lapply(1:n.perm, function(x) sample(1:N))
  S_perm_list<-lapply(sample.indexes,function(x,db) {db[x]},Sn)
  S_perm<-do.call(cbind,S_perm_list)
  
  mf<-cbind(mf,S_perm)
  
  lengths <- by(mf[,1],mf[,1],length) 
  
  
  ##################################
  # Calculate the statistic 
  ##################################
  if(k==2){
    ##################################
    # k=2
    ##################################
      #Weights for rescaling the variance
      weights<-lapply(lengths, function(x) N/x)  
      weights<-do.call(cbind,weights)
      
      
      if(test=="means"){
          mus <- apply(mf[,-1],2,function(x) by(x,mf[,1],mean) )
          sigmas <- apply(mf[,-1],2,function(x) by(x,mf[,1],var) )
        
        
          top<-sqrt(N)*diff(mus)
          bottom<-sqrt(weights%*%as.matrix(sigmas))
        
          stat<-as.numeric(top/bottom)
          parameters <- by(mf[,2],groups,mean)
          
      } else if(test=="medians"){
          mus <- apply(mf[,-1],2,function(x) by(x,mf[,1],median) )
          sigmas <- apply(mf[,-1],2,function(x) by(x,mf[,1],boot.var) )
      
          top<-sqrt(N)*diff(mus)
          bottom<-sqrt(weights%*%as.matrix(sigmas))
      
          stat<-as.numeric(top/bottom)
          parameters <- by(mf[,2],groups,median)
          
      } else if(test=="variances"){
          mus <- apply(mf[,-1],2,function(x) by(x,mf[,1],var) )
          var_sq <- apply(mf[,-1],2,function(x) by(x,mf[,1],var) )
          var_sq <-var_sq^2
          w_var_sq<-as.matrix((lengths-3)/(lengths-1),ncol=1,nrow=2)
          var_sq_w<-apply(var_sq,2,function(x) w_var_sq*x)
          f4_moment<-apply(mf[,-1],2,function(x) by(x,mf[,1],mom_4) )
          sigmas<-f4_moment-var_sq_w
      
          top<-sqrt(N)*diff(mus)
          bottom<-sqrt(weights%*%as.matrix(sigmas))
          
          
          stat<-as.numeric(top/bottom)
          parameters <- by(mf[,2],groups,var)
      }
  }else if(k>2){
    ##################################
    #  k>2
    ##################################
    weights_v<-lapply(lengths, function(x) sqrt(x))  
    weights_v<-do.call(cbind,weights_v)
    
    weights_d<-lapply(lengths, function(x) N/x)  
    weights_d<-do.call(cbind,weights_d)
    
    one<-rep(1,k)
    
    if(test=="means"){
      stat <- apply(mf[,-1],2,calc_stat_mean, groups,one,weights_v,weights_d,k)
      parameters <- by(mf[,2],groups,mean)
    }else if(test=="medians"){
      stat <- apply(mf[,-1],2,calc_stat_median, groups,one,weights_v,weights_d,k)
      parameters <- by(mf[,2],groups,median)
    }else if(test=="variances"){
      stat <- apply(mf[,-1],2,calc_stat_variance, groups,one,weights_v,weights_d,k,lengths)
      parameters <- by(mf[,2],groups,var)
    }
      
      
      
      
  }else{
    ##################################
    #  k<2
    ##################################
    warning("Need at least one group")
  }
  
  # Observed Stat
  T.obs<-stat[1]
  T.perm<-stat[-1]
  
  
  ind.rule<- mean(ifelse(T.perm>=T.obs,1,0))
  
  
  
  object_perm<-list()
  object_perm$description<-test
  object_perm$n_populations<-k
  object_perm$N<-N
  object_perm$T.obs<-T.obs
  object_perm$pvalue<-ind.rule
  object_perm$T.perm<- T.perm
  object_perm$n_perm<- n.perm
  object_perm$parameters<-c(parameters)
  object_perm$sample_sizes<-c(lengths)
  
  
  class(object_perm)<-"RPT"
  
  return(object_perm)
}


    

"boot.var" <- function(X){
  # Bootstrap method of Efron (1979). Let T(X) be the sample median
  # T(X)=X_(m) where X_(1)<=X_(2)<=,...,X_(n) is the order statistic
  # We normally assume an odd sample size n=2m-1 for convenience. We
  # will break ties by assigning m to be equal to the largest integer
  # less than (m-1)/2
  
  n=length(X)
  t=floor((n-1)/2)
  # Bootstrap Distribution
  Prob=sapply(seq(1:n), function(x) pbinom(t, size=n, prob=(x-1)/n))-sapply(seq(1:n), function(x) pbinom(t, size=n, prob=x/n))
  # Estimate of the variance of the median
  return(n*sum(((X-median(X))^2)*Prob))
}

"mom_4" <- function(X){
  #get fourth sample moment
  mean((X-mean(X))^4)
  
}



"calc_stat_mean"<-function(x,groups,one,weights_v,weights_d,k){
  mu<-as.numeric(by(x,groups,mean))
  sigmas<-as.numeric(by(x,groups,var))
  V<-as.numeric((weights_v*mu)/sigmas)
  
  d<-as.numeric(weights_d*sigmas)
  d_inv<-solve(diag(d))
  bot<-as.numeric(t(one)%*%d_inv%*%one)
  top<-sqrt(d_inv)%*%one%*%t(one)%*%sqrt(d_inv)
  top<-top/bot
  P<-diag(k)-top
  t(V)%*%P%*%V
}



"calc_stat_median"<-function(x,groups,one,weights_v,weights_d,k){
  mu<-as.numeric(by(x,groups,median))
  sigmas<-as.numeric(by(x,groups,boot.var))
  V<-as.numeric((weights_v*mu)/sigmas)
  
  d<-as.numeric(weights_d*sigmas)
  d_inv<-solve(diag(d))
  bot<-as.numeric(t(one)%*%d_inv%*%one)
  top<-sqrt(d_inv)%*%one%*%t(one)%*%sqrt(d_inv)
  top<-top/bot
  P<-diag(k)-top
  t(V)%*%P%*%V
}



"calc_stat_variance"<-function(x,groups,one,weights_v,weights_d,k,lengths){
  mu<-as.numeric(by(x,groups,var))
  var_sq <-mu^2
  w_var_sq<-as.matrix((lengths-3)/(lengths-1),ncol=1,nrow=2)
  var_sq_w<-w_var_sq*var_sq
  f4_moment<-as.numeric(by(x,groups,mom_4))
  sigmas<-as.numeric(f4_moment-var_sq_w)
  
  V<-as.numeric(weights_v*mu/sigmas)
  
  d<-as.numeric(weights_d*sigmas)
  d_inv<-solve(diag(d))
  bot<-as.numeric(t(one)%*%d_inv%*%one)
  top<-sqrt(d_inv)%*%one%*%t(one)%*%sqrt(d_inv)
  top<-top/bot
  P<-diag(k)-top
  t(V)%*%P%*%V
}
