summary.RPT<-function(object, ..., digits=max(3, getOption("digits") - 3)){

  cat("\n")
  cat("***********************************************************\n")
  cat("**   Exact and Asymptotically Robust Permutation Tests   **\n")
  cat("***********************************************************\n")

  if(object$description%in%c("means","medians","variances")) {
    cat(paste("Testing Problem: Difference of ", object$description  ,sep=""))
  } else if(object$description=="lehmann.2S.test"){
    cat("Testing Problem: Two populations differ only in location against the alternative that one population is more spread out than the other, using Lehmann (1951) two-sample U-statistic.")
  }else if(object$description=="wilcoxon.2s.test"){ 
  cat("Testing Problem: Comparison of means from two continuous distributions that satisfy a shift model assumption.")
  }else if(object$description=="hollander.2S.test"){
    cat("Testing Problem: Comparison of means from two continuous distributions that satisfy a shift model assumption, using Hollander (1967) two-sample U-statistic.")
  }
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
    }else if(object$description=="lehmann.2S.test"){
      description_2<- "Pr(|Y-Y'|>|X-X'|)=1/2" 
    }else if(object$description=="wilcoxon.2s.test"){
      if(object$wilcoxon.type=="continuity") description_2<-"Pr(X<=Y)=1/2"
      if(object$wilcoxon.type=="discontinuity") description_2<-"Pr(X<=Y)=Pr(Y<=X)"
    }else if(object$description=="hollander.2S.test"){
      description_2<-" Pr(X+X'<Y+Y')=1/2"
    }
  cat("* --------------------------------------------------------*\n")
  cat(paste("H0: ", description_2 ,sep=""))
  cat("\n")
  cat("* --------------------------------------------------------*\n")
  cat("Estimates:\n")
  cat("\n")
  if(object$description %in% c("means","medians","variances")) {
    z<-cbind(object$parameters,object$sample_sizes)
    colnames(z)<-c("Parameter", "Sample Size")
    print(z)
    cat("--------\n")
  }
  cat(paste("Test Statistic: ", round(object$T.obs,3) ,sep=""))
  cat("\n")
  cat(paste("P-Value: ", round(object$pvalue,3) ,sep=""))
  cat("\n")

}


