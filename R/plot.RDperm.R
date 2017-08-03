#' @title Plot RDperm
#' @description Plots a histogram and empirical cdf
#'
#'
#' @method plot RDperm
#'
#'
#' @param x Object of class "RDperm"
#' @param w Character. Name of variable to be plotted
#' @param plot.class Character. Can be: "both" for a histogram and cdf plot, "hist" for a histogram or "cdf" for only the cdf plot
#' @param ... Additional ggplot2 controls
#'
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Canay, I and Kamat V. (2016) ``Approximate Permutation Tests and Induced Order Statistics in the Regression Discontinuity Design''
#' @keywords permutation rdperm
#' @importFrom stats ecdf
#' @importFrom graphics plot
#' @import ggplot2
#' @import gridExtra
#' @examples
#' permtest<-RDperm(W=c("demshareprev"),z="difdemshare",data=lee2008)
#' plot(permtest,w="demshareprev")
#' @export

plot.RDperm <- function(x,w, plot.class="both",...){
  if(class(x)!="RDperm")(print("Element has to be of class RDperm"))
  data<-x$data
  cutoff<-x$cutoff
  z<-x$rv
  q<-x$results[rownames(x$results)==w,3]
  W_z<-base::subset(data, select=c(w,z))
  colnames(W_z)[colnames(W_z)==z]<-"z"
  N<-dim(data)[1]

  W_left <- W_z[W_z$z<cutoff,]
  n_left <- length(W_left$z)
  W_right <- W_z[W_z$z>=cutoff,]
  n_right <- length(W_right$z)

  # Induced order of W obs
  W_left <- W_left[order(W_left$z),]
  W_right <- W_right[order(W_right$z),]

  W_left_q<-base::subset(W_left[(n_left-q+1):n_left,], select=c(w))
  Z_left<-base::subset(W_left[(n_left-q+1):n_left,], select=c(z))
  W_right_q<-base::subset(W_right[1:q,], select=c(w))
  Z_right <-base::subset(W_right[1:q,], select=c(z))


  # Create the data frames.
  # I: for the conditional mean
  dat <- data.frame(rbind(W_left_q,W_right_q),rbind(Z_left,Z_right))
  colnames(dat) <- c(w,"Z")

  # II: For the ECDF and Histograms
  left<-data.frame(variable=rep("Left of threshold",q), value=W_left_q)
  right<-data.frame(variable=rep("Right of threshold",q), value=W_right_q)
  gdat <- data.frame(rbind(left,right))
  colnames(gdat)<-c("variable","value")
  # Create the histograms using ggplot.
  hist <- ggplot2::ggplot(gdat, aes_string(x="value", fill="variable")) +
                    geom_histogram(alpha=0.2, position="identity", bins=30) +
                    theme_bw() +
                    theme(legend.justification=c(0,0),
                          legend.position=c(0.01,0.8),
                          panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
                    labs(fill="",y="Density",x="") +
                    scale_color_manual(labels = c("Left of threshold", "Right of threshold"),
                                       values = c("#FFFF00", "#9933FF"))

  # Create the ECDFs using ggplot.
  cdf <- ggplot2::ggplot(gdat, aes_string(x="value")) +
          stat_ecdf(aes_string(colour="variable")) +
          theme_bw() +
          theme(legend.justification=c(0,0),
                legend.position=c(0.01,0.8),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          labs(fill="",y="Emp. CDF",x="",color="") +
          scale_color_manual(labels = c("Left of threshold", "Right of threshold"),
                            values = c("#FF9900", "#9933FF"))



  if(plot.class=="hist"){
    return(hist)
  }else if(plot.class=="cdf"){
    return(cdf)
  }else{
    return(gridExtra::grid.arrange(hist, cdf, ncol=2))}

}
