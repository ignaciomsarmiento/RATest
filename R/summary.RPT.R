
#' @export



summary.RPT<-function(object, ..., digits=max(3, getOption("digits") - 3)){

  cat("\n")
  cat("***********************************************************\n")
  cat("**   Exact and Asymptotically Robust Permutation Tests   **\n")
  cat("***********************************************************\n")

  cat(paste("Testing Problem: Difference of ", object$description  ,sep=""))
  cat("\n")
  cat(paste("Number of Populations: ", object$n_populations ,sep=""))
  cat("\n")
  cat(paste("Total Number of Observations (Pooled Sample): ", object$N ,sep=""))
  cat("\n")
  cat(paste("Number of Permutations: ",object$n_perm,sep=""))
  cat("\n\n")
  if(object$description=="means"){
    description_2<-"Population means are equal"
    }else if(object$description=="medians"){
    description_2<-"Population medians are equal"
    }else if(object$description=="variances"){
    description_2<-"Population variances are equal"
    }
  cat("* --------------------------------------------------------*\n")
  cat(paste("H0: ", description_2 ,sep=""))
  cat("\n")
  cat("* --------------------------------------------------------*\n")
  cat("Estimates:\n")
  cat("\n")
  z<-cbind(object$parameters,object$sample_sizes)
  colnames(z)<-c("Parameter", "Sample Size")
  print(z)
  cat("--------\n")
  
  cat(paste("Test Statistic: ", round(object$T.obs,3) ,sep=""))
  cat("\n")
  cat(paste("P-Value: ", object$pvalue ,sep=""))
  cat("\n")

}




#
#
#   * THIS IS THE DICTIONARY OF DESCRIPTIONS FOR EVERY CASE CONSIDERED ABOVE
#
#   * When it's equality of means (any number of samples)
#
#   DESCRIPTION_1 Comparison of means from multiple populations.
#   DESCRIPTION_2 Population means are equal
#
#   * When it's equality of variances (any number of samples)
#
#   DESCRIPTION_1 Comparison of variances from multiple populations.
#   DESCRIPTION_2 Population variances are equal.
#
#
#   * When it's equality of medians (any number of samples)
#
#   DESCRIPTION_1 Comparison of medians from multiple populations.
#   DESCRIPTION_2 Population medians are equal.
#
#
#   * When it's Lehmann's U statistic
#
#   DESCRIPTION_1 Two populations differ only in location against the alternative that one population is more spread out than the other, using Lehmann (1951) two-sample U-statistic.
#   DESCRIPTION_2 Pr(|Y-Y'|>|X-X'|)=1/2
#
#
#   * When it's Wilcoxon two-sample Statistic
#
#   DESCRIPTION_1 Comparison of means from two continuous distributions that satisfy a shif model assumption.
#   DESCRIPTION_2 Pr(X<=Y)=1/2
#
#
#   * When it's Wilcoxon Statistic without continuity assumption
#
#   DESCRIPTION_1 Comparison of means from two distributions that satisfy a shif model assumption.
#   DESCRIPTION_2 Pr(X<=Y)=Pr(Y<=X)
#
#
#   * When it's Hollander's U statistic
#
#   DESCRIPTION_1 Comparison of means from two continuous distributions that satisfy a shif model assumption, using Hollander (1967) two-sample U-statistic.
#   DESCRIPTION_2 Pr(X+X'<Y+Y')=1/2
