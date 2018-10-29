#' @title Robust Permutation Test
#'
#' @description 
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the groups on the right.
#' @param data a data.frame in which to interpret the variables named in the formula. If this is missing, then the variables in the formula should be on the search list. 
#' @param test test to be perfomed. Multiple options are available, depending on the nature of the testing problem. In general, we have two types of problem. First, when the researcher is interested in comparing parameters. In this case, "means" will perform a Difference of Means, "medians" a Difference of Medians, "variances" a Difference of Variances. This case allows for 2 or more population comparisons. For the test of difference of medians the Efron (1992) bootstrap estimator is used to estimate the variances (for further details, see Chung and Romano (2013)). Second, when the parameter of interest is a function of the joint distribution. In this case, "lehmann.2S.test" will perform Lehmann (1951) two-sample U statistics, "wilcoxon.2s.test" the two-sample Wilcoxon test (with or without continuity assumption), and "hollander.2S.test" Hollander (1967) two sample U statistics. In this case, only 2 sample comparisons are permitted. 
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. See remark 3.2 in Canay and Kamat (2017). The default is n.perm=499.	
#' @param na.action a function to filter missing data. This is applied to the model.frame . The default is na.omit, which deletes observations that contain one or more missing values.

#' @return An object of class "RPT" is a list containing at least the following components:


#'  \item{description}{Type of test, can be Difference of Means, Medians, or Variances.}
#'  \item{n_populations}{Number of grups.}
#'  \item{N}{Sample Size.}
#'  \item{T.obs}{Observed test statistic.}
#'  \item{pvalue}{P-value.}
#'  \item{T.perm}{Vector. Test statistics from the permutations.}
#'  \item{n_perm}{Number of permutations.}
#'  \item{parameters}{Estimated parameters.}
#'  \item{sample_sizes}{Groups lengths.}

#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references 
#' Chung, E. and Romano, J. P. (2013). Exact and asymptotically robust permutation tests. The Annals of Statistics, 41(2):484–507.
#' Chung, E. and Romano, J. P. (2016). Asymptotically valid and exact permutation tests based on two-sample u-statistics. Journal of Statistical Planning and Inference, 168:97–105.
#' Devroye, L. P. and Wagner, T. J. (1980). The strong uniform consistency of kernel density estimates. In Multivariate Analysis V: Proceedings of the fifth International Symposium on Multivariate Analysis, volume 5, pages 59–77.
#' Efron, B. (1992). Bootstrap methods: another look at the jackknife. In Breakthroughs in statistics, pages 569–593. Springer.
#' Hall, P., DiCiccio, T. J., and Romano, J. P. (1989). On smoothing and the bootstrap. The Annals of Statistics, pages 692–704.
#' Hollander, M. (1967). Asymptotic efficiency of two nonparametric competitors of wilcoxon’s two sample test. Journal of the American Statistical Association, 62(319):939–949.
#' Lehmann, E. L. (1951). Consistency and unbiasedness of certain nonparametric tests. The Annals of Mathematical Statistics, pages 165–179.
#' @keywords robust permutation test rpt
#' @import quantreg
#' @importFrom stats cor var runif
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



RPT<-function(formula,data,test="means",n.perm=499, na.action, wilcoxon.option="continuity"){
  if(!(test%in%c("means","medians","variances","lehmann.2S.test","wilcoxon.2s.test","hollander.2S.test"))){
    print("Must specify test to perform. Options include 'means', 'medians', 'variances', 'Lehmann's two-sample U statistics','Two-sample Wilcoxon test' (with or without continuity assumption), or 'Hollander's two sample U statistics' ")
    stop()
  }
  #Put the formula into a model.frame
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data",  "na.action"), names(mf), 0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  
  #puts the categories into factors 
  mf[,2] <- as.factor(mf[,2])
    
  N<-dim(mf)[1] #sample size
  
  k<-length(levels(mf[,2])) #number of samples 
  
  names<-as.character(unique(mf[,2])) #names of the categories
  
  #Permutations
  Sn<-mf[,1] #Sn are the stacked values of the outcome
  groups<-mf[,2] #are the groups
  mf<-mf[,c(2,1)]#change the order 
  
  # Generate random permutations
  sample.indexes = lapply(1:n.perm, function(x) sample(1:N))
  S_perm_list<-lapply(sample.indexes,function(x,db) {db[x]},Sn)
  S_perm<-do.call(cbind,S_perm_list)
  
  mf<-cbind(mf,S_perm) #first column is the group indicator, second original/observed data, rest are permutations
  
  lengths <- by(mf[,1],mf[,1],length)  #observations per group
  
  
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
          
      } else if(test=="lehmann.2S.test"){
          m<-lengths[[1]]
          n<-lengths[[2]]
          stat <- apply(mf[,-1],2, function(x) lehmann.2S.test(x,m,n))
          
      } else if(test=="wilcoxon.2s.test" ){
        m<-lengths[[1]]
        n<-lengths[[2]]
        if(wilcoxon.option=="discontinuity"){
          stat <- apply(mf[,-1],2, function(x) wilcoxon.2s.test(x,m,n,type="discontinuity"))
        }else if(wilcoxon.option=="continuity"){
          stat <- apply(mf[,-1],2, function(x) wilcoxon.2s.test(x,m,n,type="continuity"))
        } 

      } else if(test=="hollander.2S.test"){
        m<-lengths[[1]]
        n<-lengths[[2]]
        stat <- apply(mf[,-1],2, function(x) hollander.2S.test(x,m,n))
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
  
  
  T.obs<-stat[1] # Observed Stat
  T.perm<-stat[-1] #Permutations
  
  #Indicator rule
  ind.rule<- mean(ifelse(T.perm>=T.obs,1,0))
  
  
  
  object_perm<-list() #Generates an empty list to collect all the required info for summary
  object_perm$description<-test
  if(test=="wilcoxon.2s.test") object_perm$wilcoxon.type<-wilcoxon.option
  object_perm$n_populations<-k
  object_perm$N<-N
  object_perm$T.obs<-T.obs
  object_perm$pvalue<-ind.rule
  object_perm$T.perm<- T.perm
  object_perm$n_perm<- n.perm
  if(test %in% c("means","medians","variances")) object_perm$parameters<-c(parameters)
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

"lehmann.2S.test" <- function(perX,m,n){
  N<-m+n
  X <- perX[1:m]
  Y <- perX[(m+1):N]
  if (anyNA(c(X,Y))) 
    stop("NAs in first or second argument")
  if (!is.numeric(c(X,Y))) 
    stop("Arguments must be numeric")
  else
    m = length(X)
    n = length(Y)
    W = lapply(list(X,Y), function(x) abs(as.vector(t(outer(x,x, FUN="-")))))
    Z = outer(W[[1]],W[[2]], FUN=">")
  
    V = 4*(  (m-1)^{-1}*sum (unlist(lapply(1:(m-1),function(i) quadratic.sum(Z,m,n,i)))) + 
             (m/n)*(n-1)^{-1}*sum (unlist(lapply(1:(n-1),function(i) quadratic.sum(t(Z),n,m,i))))   )
  
    return( sum(Z-0.5)/(m*n)^2*sqrt(V))
}




"wilcoxon.2s.test" <- function(perX,m,n, type=c("continuity","discontinuity")){
  N<-m+n
  X <- perX[1:m]
  Y <- perX[(m+1):N]
  if (anyNA(c(X,Y))) stop("NAs in first or second argument")
  if (!is.numeric(c(X,Y)))  stop("Arguments must be numeric")
  
  
  type = match.arg( type )
  
  # generate Y0
  if(type == "continuity"){
    Z = outer(X,Y, FUN="<=")
    V = (m*(m-1))^{-1}*sum(unlist(lapply(1:m, function(l) (mean(Z[l,])- m^{-1}*sum(unlist(lapply(1:m, function(i) mean(Z[i,])))))^2))) +
      (n*(n-1))^{-1}*sum(unlist(lapply(1:n, function(l) (mean(Z[,l])- n^{-1}*sum(unlist(lapply(1:n, function(i) mean(Z[,i])))))^2)))
    Tmn = (sum(Z)-(m*n)/2)/(m*n*sqrt(V))
  }
  else if(type == "discontinuity"){
    Z.l = outer(X,Y, FUN="<")
    Z.e = outer(X,Y, FUN="==")
    
    V = (m*(m-1))^{-1}*sum(unlist(lapply(1:m, function(l) (mean(Z.l[l,]+0.5*Z.e[l,])- m^{-1}*sum(unlist(lapply(1:m, function(i) mean(Z.l[i,]+0.5*Z.e[l,])))))^2))) +
      (n*(n-1))^{-1}*sum(unlist(lapply(1:n, function(l) (mean(Z.l[,l]+0.5*Z.e[,l])- n^{-1}*sum(unlist(lapply(1:n, function(i) mean(Z.l[,i]+0.5*Z.e[,l])))))^2)))
    Tmn = (sum(Z.l)+.5*sum(Z.e)-(m*n)/2)/(m*n*sqrt(V))
    
  } else {
    stop( paste( "Unrecognized type '", type, "'", sep="" ) )
  }
  
  return( Tmn )
}



"hollander.2S.test" <- function(perX,m,n){
  N<-m+n
  X <- perX[1:m]
  Y <- perX[(m+1):N]
  if (anyNA(c(X,Y))) 
    stop("NAs in first or second argument")
  if (!is.numeric(c(X,Y))) 
    stop("Arguments must be numeric")
  else
    m = length(X)
  n = length(Y)
  W = lapply(list(X,Y), function(x) as.vector(t(outer(x,x, FUN="+"))))
  Z = outer(W[[1]],W[[2]], FUN=">")
  
  V = 4*(  (m-1)^{-1}*sum (unlist(lapply(1:(m-1),function(i) quadratic.sum(Z,m,n,i)))) + 
             (m/n)*(n-1)^{-1}*sum (unlist(lapply(1:(n-1),function(i) quadratic.sum(t(Z),n,m,i))))   )
  
  return( sum(Z-0.5)/(m*n)^2*sqrt(V))
}


quadratic.sum <- function(Z,m,n,i){
  A = (csi(Z,m,n,i)-(m-1)^{-1}*sum(unlist(lapply(1:(m-1), function(i) csi(Z,m,n,i)))))^2
  return(A)
}

csi <- function(Z,m,n,i){
  C = lapply((i+1):m, function(j) lapply(1:(n-1), function(k) lapply((k+1):n, function(l) Z[((i-1)*m+j):m,((k-1)*n+l):(n+(k-1)*n)])))
  return( 2/((m-i)*n*(n-1))*sum(unlist(C)) )
}

