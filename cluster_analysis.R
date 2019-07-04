


#' Cluster level estimate of the treatment effect as described by Hughes in
#' "Robust inference for the stepped wedge design"
#' Computes the treatment effect in terms of a risk difference.
#'
#' @export
#' @param d data from stepped wedge design
#' @return point estimate of the treatment effect
delta <- function(d){
  
  d1 <- d[, .(Y = mean(y), x = mean(tx)), keyby = .(clust_id, t)]
  d2 <- d[, .(xhat = mean(tx)), keyby = .(t)]
  d1 <- merge(d1, d2, by = c("t"))
  
  d1[, numer:= Y*(x - xhat)]
  d1[, denom:= xhat*(1-xhat)]
  
  sum(d1$numer)/(length(unique(d1$clust_id))*sum(d1$denom))
}
