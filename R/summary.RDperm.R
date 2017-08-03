
#' @export



summary.RDperm<-function(object, ..., digits=max(3, getOption("digits") - 3)){

  cat("\n")
  cat("**********************************************************\n")
  cat("**       RD Distribution Test using permutations        **\n")
  cat("**********************************************************\n")

  cat(paste("Running Variable: ",object$rv,sep=""))
  cat("\n")
  cat(paste("Cutoff: ",object$cutoff,sep=""))
  cat("\n")
  cat(paste("q: ",object$q_type,sep=""))
  cat("\n")
  cat(paste("Test Statistic: ",object$test.statistic,sep=""))
  cat("\n")
  cat(paste("Number of Permutations: ",object$n_perm,sep=""))
  cat("\n")
  cat(paste("Number of Obs: ", dim(object$data)[1],sep=""))

  cat("\n\n")
  cat("**********************************************************\n")
  cat("H0: 'Continuity of the baseline covariates at the cutoff' \n")
  cat("**********************************************************\n\n")
  cat("Estimates:\n")
  n<-dim(object$results)[1]
  stars<-vector(length=n)
  object$results[,1:2]<-round(object$results[,1:2],2)
  for(i in 1:n) {
    stars[i]<- if(object$results[i,2]<0.01) "***" else if(object$results[i,2]<0.05) "**" else if(object$results[i,2]<0.01) "*" else " "
  }


  print.default(cbind(object$results," "=stars),quote=FALSE,print.gap=2,right=FALSE)

  cat("---\n")
  cat("Signif. codes:   0.01 '***' 0.05 '**' 0.1 '*' \n\n")
  results<-list(results=object$results)
  class(results)<-"summary.RDperm"
  return(invisible(results))
}
