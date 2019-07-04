library(data.table)
library(splitstackshape)
library(simstudy)
library(ggplot2)

#' Computes the cluster level variance for cohort SW CRT based on specified
#' between subject correlation and between cluster correlation
#' See file:///C:/Users/mjones/Documents/examples/000_simulation/icc/001_icc_3_level_binary.html
#'
#' @export
#' @param rho2 is the correlation between participants within a cluster
#' @param rho3 is the proportion of total variance which due to between cluster variance
#' @return var_k
#' @example 
#' rho2 <- 0.3
#' rho3 <- 0.05
#' var_clust(rho2, rho3)
#' [1] 0.1973921
var_clust <- function(rho2, rho3){
  ((pi^2)/3) * (rho3*rho2/(rho2-rho3))
}


#' Computes the subject level variance component for cohort SW CRT based on 
#' specified between subject correlation and between cluster correlation.
#' See file:///C:/Users/mjones/Documents/examples/000_simulation/icc/001_icc_3_level_binary.html
#'
#' @export
#' @param rho2 is the correlation between participants within a cluster
#' @param rho3 is the proportion of total variance which due to between cluster variance
#' @return var_subj
#' @example 
#' rho2 <- 0.3
#' rho3 <- 0.05
#' var_subj(rho2, rho3)
#' [1] 0.4605815
var_subj <- function(rho2, rho3){
  ((pi^2)/3) * (rho3*rho2/(rho2-rho3)) * ((1/rho2)-1)
}











#' Computes the 'a' parameter for the beta distribution to be used to generate
#' the baseline cluster proportion for a cross sectional SW CRT.
#' See file:///C:/Users/mjones/Documents/examples/000_simulation/icc/001_icc_3_level_binary.html
#'
#' @export
#' @param rho desired ICC
#' @param mu mean 
#' @return a
#' @example 
#' a(rho = 0.05, mu = 0.6)
#' [1] 11.4
a <- function(rho, mu)  (1-rho)*mu/rho

#' Computes the 'b' parameter for the beta distribution to be used to generate
#' the baseline cluster proportion for a cross sectional SW CRT.
#' See file:///C:/Users/mjones/Documents/examples/000_simulation/icc/001_icc_3_level_binary.html
#'
#' @export
#' @param rho desired ICC
#' @param mu mean 
#' @return a
#' @example 
#' b(rho = 0.05, mu = 0.6)
#' [1] 7.6
b <- function(rho, mu)  (1-rho)*(1-mu)/rho

#' Converts probability to odd scale.
#'
#' @param p probability value
#' @return odd
#' @examples 
#' odd(0.4)
#' [1] 0.6666667
#' 
#' odd(0.5)
#' [1] 1
odd <- function(p) p/(1-p)


#' Computes logit (log-odds)
#'
#' @param p probability value
#' @return log-odds
#' @examples 
#' logit(0.4)
#' [1] -0.4054651
#' 
#' logit(0.5)
#' [1] 0
#' 
#' logit(0.6)
#' [1] 0.4054651
logit <- function(p) log(odd(p))


#' Computes inverse logit giving probability. 
#'
#' @param p probability value
#' @return log-odds
#' @examples 
#' inv_logit(0)
#' [1] 0.5
#' 
#' inv_logit(0.4054651)
#' [1] 0.6
inv_logit <- function(x) exp(x)/(1+exp(x))


#' Generates SW CRT dataset containing cluster ID, timepoints, transition
#' start time, current treatment status, subject ID etc.
#' Each cluster has a specific baseline probability of event which is 
#' used in computing the subject level probability of event that is ultimately
#' used to generate a bernoulli draw for to given an event status (0/1).
#' The linear predictor includes terms for the baseline probability of an event
#' the treatment status and a (linear) temporal trend term.
#'
#' @export
#' @param ldes_spec Design specification
#' @return long data.frame of subject level responses
dgp_cross_sec <- function(ldes_spec){
  
  # make a big dataset and then sample from it
  
  d <- data.table::CJ(clust_id = 1:ldes_spec$n_clusters,
                      t = ldes_spec$t_obs)
  
  gen_n_sub <- function(idx){
    
    nmin <- ldes_spec$nmin_clust[idx]
    nmax <- ldes_spec$nmax_clust[idx]
    
    sample(nmin:nmax, 
           size = length(ldes_spec$t_obs), 
           replace = T)
    
  }
  
  d$n_subj <- unlist(lapply(1:ldes_spec$n_clusters, gen_n_sub))
  
  # regenerate the data.table now with the correct number of rows for subjects
  d <- splitstackshape::expandRows(d, "n_subj", count.is.col = T, drop = F)
  
  
  for(i in 1:ldes_spec$n_clusters){
    
  }
  
  
  
  # n for each timepoint for each cluster
  n_subj <- sample(x = ldes_spec$nmin_clust:ldes_spec$nmax_clust, 
                   size = ldes_spec$n_clusters*length(ldes_spec$t_obs), 
                   replace = T)
  
  gen_sub <- function(idx){
    
    data.table::CJ(clust_id = 1:ldes_spec$n_clusters,
                   t = ldes_spec$t_obs,
                   subj_id = )
    
    
  }
  
  lapply(1:length(n_subj), gen_sub)
  
  d <- data.table::CJ(clust_id = 1:ldes_spec$n_clusters,
                      t = ldes_spec$t_obs,
                      subj_id = )
  
  # data for cluster, timepoint and individual
  d <- data.frame(expand.grid(clust_id = 1:n_clust,
                              t = 0:n_time,
                              sub_id = 1:n_sub),
                  stringsAsFactors = F)
  
  d <- merge(d, design_mat, by = "clust_id")
  d$tx_active <- as.numeric(d$t >= d$tx_start)
  # just to get unique sub ids
  d$sub_id_full <- paste("k", d$clust_id, "t", d$t, "i", d$sub_id, sep = "")
  d <- d[order(d$clust_id, d$t, d$sub_id),]
  # beta binomial
  d$p0 <- rep(rbeta(n_clust, shape1 = a(rho, mu), shape2 = b(rho, mu)),
              each = (n_time +1 ) * n_sub)
  
  # linear time trends, nothing fancy:
  # intercept
  # tx 
  d$eta <- log(odd(d$p0)) + beta1 * d$tx_active + beta2 * d$t
  d$p <- inv_logit(d$eta)
  d$y <- rbinom(nrow(d), 1, prob = d$p)
  
  # initially set the sequence order to the cluster id
  d$seq_num <- d$clust_id 
  d
}




dgp_closed_cohort <- function(ldes_spec){
  
  # variance for re
  var_c <- var_clust(ldes_spec$rho2, ldes_spec$rho3)
  var_i <- var_subj(ldes_spec$rho2, ldes_spec$rho3)
  
  # observations occur at following months
  dlu1 <- data.table(t = ldes_spec$t_obs,
                     mnth = ldes_spec$months)
  
  # tx start times
  dlu2 <- data.table(clust_id = 1:ldes_spec$n_clusters,
                     tx_start = ldes_spec$tx_start)
  
  # clust intercepts
  dlu3 <- data.table(clust_id = 1:ldes_spec$n_clusters,
                     clust_int = rnorm(ldes_spec$n_clusters, 0, var_c))
  
  # subject level intercepts (subj have repeat measures)
  dlu4 <- data.table::CJ(clust_id = 1:ldes_spec$n_clusters,
                         subj_id = 1:ldes_spec$n_subj)
  dlu4[,subj_int := rnorm(nrow(dlu4), 0, var_i)]
  
  
  
  # basic data.table 
  d <- data.table::CJ(clust_id = 1:ldes_spec$n_clusters,
                      subj_id = 1:ldes_spec$n_subj,
                      t = ldes_spec$t_obs)
  
  d[,subj_id_full:= paste0(clust_id, ":", subj_id)]
  
  # merge in month at which visit occurs
  d <- merge(d, dlu1, by = "t")
  # merge in tx start time
  d <- merge(d, dlu2, by = "clust_id")
  # add tx indicator
  d[,tx:=as.numeric(t>=tx_start)]
  # merge in cluster level intercept
  d <- merge(d, dlu3, by = "clust_id")
  # merge in subject level intercept for repeat measures
  d <- merge(d, dlu4, by = c("clust_id", "subj_id"))
  # add a cluster sub concatenation
  d[, cs_id := paste0(clust_id, ":",subj_id)]
  
  # fixef
  log_odds_0 <- log(odd(ldes_spec$p_baseline))
  log_ortrt <- log(odd(ldes_spec$p_tx)/odd(ldes_spec$p_baseline))
  log_ortime <- log(odd(ldes_spec$p_t_effects)/odd(ldes_spec$p_baseline))
  
  d$eta <- log_odds_0 + d$clust_int + d$subj_int + log_ortrt * d$tx +
    log_ortime[1] * as.numeric(d$t == 1) +
    log_ortime[2] * as.numeric(d$t == 2) +
    log_ortime[3] * as.numeric(d$t == 3) +
    log_ortime[4] * as.numeric(d$t == 4) +
    log_ortime[5] * as.numeric(d$t == 5) +
    log_ortime[6] * as.numeric(d$t == 6) +
    log_ortime[7] * as.numeric(d$t == 7) +
    log_ortime[8] * as.numeric(d$t == 8)
  
  d$p <- inv_logit(d$eta)
  d$y <- rbinom(nrow(d), 1, d$p)
  d
  
}