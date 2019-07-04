library(lme4)
library(gee)
library(brms)
library(rstan)
library(MCMCvis)

rstan::lookup("inv_logit")
rstan::rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = parallel::detectCores())


source("dgp.R")
source("cluster_analysis.R")
# see file:///C:/Users/mjones/Documents/examples/000_simulation/icc/001_icc_3_level_binary.html
# generate data for a cross-sectional SW design with a known ICC you can generate the 
# baseline cluster proportions from a beta distribution, p0k∼Beta(a,b).


# Design specification for the sw crt
# For closed cohort we observe the same people from a given intercept at
# each timepoint - i.e. we have repeat measures.
# n_clusters number of clusters
# n_subj number of subjects observed at each cluster at each timepoint 
# months the months that observations take place
# t_obs vector of observation visits
# t_obs_baseline observations that are baseline measures 
# t_obs_fu observations that are follow up (post intervention)
# tx_start start time for tx in each cluster
# p_baseline baseline prob of event
# p_tx prob of event under active treatment 
# p_t_effects temporal prob of event for ctl grp (from first obs i.e. month 3)
# rho2 between person correlation (within cluster)
# rho3 between cluster correlation
ldes_spec <- list(n_clusters = 4,
                  n_subj = 50, 
                  months = c(0,3,6,12,15,18,24,27,30),
                  t_obs = c(0,1,2,3,4,5,6,7,8),
                  t_obs_baseline = c(0,1),
                  t_obs_fu = 8,
                  tx_start = c(2,3,5,6),
                  p_baseline = 0.45,
                  p_tx = 0.2,
                  p_t_effects = 0.45 - c(3, 6, 12, 15, 18, 24, 27, 30)*0.002,
                  rho2 = 0.2,  
                  rho3 = 0.01)

# ldes_spec <- list(n_clusters = 40,
#                   n_subj = 50,
#                   months = c(0,3,6,12,15,18,24,27,30),
#                   t_obs = c(0,1,2,3,4,5,6,7,8),
#                   t_obs_baseline = c(0,1),
#                   t_obs_fu = 8,
#                   tx_start = c(2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,
#                                5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6),
#                   p_baseline = 0.45,
#                   p_tx = 0.2,
#                   p_t_effects = 0.45 - c(3, 6, 12, 15, 18, 24, 27, 30)*0.002,
#                   rho2 = 0.2,
#                   rho3 = 0.01)

d <- dgp_closed_cohort(ldes_spec)
lm1 <- glmer(y ~ 1 + tx + mnth + (1|clust_id) + (1|subj_id_full) + factor(t) , 
                     data = d,
                     family = binomial,
                     control = glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=100000),
                                            check.conv.singular="warning"))

ld <- brms::make_standata(y ~ 0 + tx + mnth + (1|clust_id) + (1|subj_id_full) + factor(t), 
                          data = d, 
                          family = bernoulli())

# factor for time was introduced erroneously
ld$X <- ld$X[,-3]
ld$K <- ncol(ld$X)

lm1 <- glm(y ~ 1 + tx + mnth + as.factor(t), data = d, family = binomial)
summary(lm1)


myf <- file.path("model001.stan")
mod <- rstan::stan_model(myf ,verbose = F)

# for diagnostics I initially had chains = 4, but changed this to 1 after i was convinced 
# that everything looked ok.
set.seed(315430)
fit1 <- rstan::sampling(object  = mod,
                       data    = ld,
                       chains  = 6,
                       thin    = 2,
                       iter    = 3000,
                       refresh = 1)
summary(fit1, pars = c("b_Intercept", "b", "sd_1", "sd_2"),
        prob = c(0.025, 0.975))$summary

posterior <- as.array(fit1)
lp_p <- bayesplot::log_posterior(fit1)
np_p <- bayesplot::nuts_params(fit1)

bayesplot::mcmc_trace(posterior, 
                      regex_pars = c( "b"), 
                      np = np_p) + 
  xlab("Post-warmup iteration")

bayesplot::mcmc_acf(posterior, 
                    pars = c( "b[1]"),  lags = 10)

MCMCvis::MCMCsummary(fit1, round = 3 , 
            params = c("b_Intercept","b","sd_1","sd_2"), 
            Rhat = F, 
            n.eff = TRUE)


nsim <- 1000
m <- matrix(0, nrow = nsim, ncol = 1)
for(i in 1:nsim){
  d <- dgp_closed_cohort(ldes_spec)
  
  m[i, ] <- delta(d)
}
hist(m[,1])

par(mfrow = c(1, 2))
hist(m[,1])
abline(v = 0.45, col = "red")
abline(v = mean(m[,1]), col = "dodgerblue")
hist(m[,2])
abline(v = 0.2, col = "red")
abline(v = mean(m[,2]), col = "dodgerblue")
par(mfrow = c(1, 1))

d <- dgp_closed_cohort(ldes_spec)
summary(d)

lm1 <- glm(y ~ tx + factor(t) , data = d, family = "binomial")
summary(lm1)

dnew <- d[1:2,]
dnew[,t:= rep(0, 2)]
dnew[,tx:= 0:1]
predict(lm1, type = "response", newdata = dnew)

d <- dgp_closed_cohort(ldes_spec)
ggplot(d, aes(x = mnth, y = p, 
               group = paste0(clust_id, ":", subj_id),
               color = paste0(clust_id)))+
  geom_line()+
  scale_x_continuous("Months from start")+
  scale_y_continuous("Probability of evnt")

dt <- d[, .(phat = mean(y)), keyby = .(clust_id, mnth, t, tx)]
ggplot(dt, aes(x = mnth, y = phat, 
              group = paste0(clust_id),
              color = paste0(clust_id)))+
  geom_line()+
  scale_x_continuous("Months from start")+
  scale_y_continuous("Probability of evnt")




# nsim <- 1000
# m <- matrix(0, nrow = nsim, ncol = 2)
# dnew <- d[1:2,]
# dnew[,t:= rep(0, 2)]
# dnew[,tx:= 0:1]
# 
# for(i in 1:nsim){
#   d <- dgp_closed_cohort(ldes_spec)
#   # summary(d)
#   
#   lm1 <- glm(y ~ tx + factor(t) + factor(clust_id) , data = d, family = "binomial")
#   # summary(lm1)
#   
#   p <- predict(lm1, type = "response", newdata = dnew)
#   m[i, ] <- p
# }
# ldes_spec <- list(n_clusters = 10,
#                   n_subj = 50, 
#                   months = c(0,3,6,12,15,18,24,27,30),
#                   t_obs = c(0,1,2,3,4,5,6,7,8),
#                   t_obs_baseline = c(0,1),
#                   t_obs_fu = 8,
#                   tx_start = c(2,2,2,2,2,5,5,5,5,5),
#                   p_baseline = 0.45,
#                   p_tx = 0.2,
#                   p_t_effects = 0.45 - c(3, 6, 12, 15, 18, 24, 27, 30)*0.002,
#                   rho2 = 0.2,  
#                   rho3 = 0.01) 
