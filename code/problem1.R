setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- p1_prelim
# Load data
library(ggplot2)
library(tibble)
source("../data/probAhelp.R")
source("../data/probAdata.R")
set.seed(123)

# Least sum of absolute residuals and least squares estimator
beta = ARp.beta.est(data3A$x, 2)
beta_la = beta$LA
beta_ls = beta$LS

# Centered residuals
epsilon_la = ARp.resid(data3A$x, beta_la)
epsilon_ls = ARp.resid(data3A$x, beta_ls)

## ---- p11_sim
# AR residual resampling method
AR2.rrbootstrap = function(x, beta, e){
  i = sample(1:99, 1)
  x0 = x[i:(i + 1)]
  return(ARp.filter(x0, beta, sample(e, replace=T)))
}

# Bootstrap samples for LA and LS estimators
B = 2000
m = length(data3A$x)
x_hat_la = t(replicate(
  B, AR2.rrbootstrap(data3A$x, beta_la, epsilon_la)))
x_hat_ls = t(replicate(
  B, AR2.rrbootstrap(data3A$x, beta_ls, epsilon_ls)))

# Plot some of the samples and original data
data = tibble(x=data3A$x, t=1:100)
p_x = ggplot(data, aes(t, x, col="Data")) + geom_line()
nums = c(300, 600, 1000)
p_time_series = p_x
for(i in 1:length(nums)){
  p_time_series = p_time_series +
    geom_line(
      data=tibble(x=x_hat_la[nums[i],], t=1:100),
      (aes(col="LA"))) +
    geom_line(
      data=tibble(x=x_hat_ls[nums[i],], t=1:100),
      (aes(col="LS"))) +
    labs(col="")
}
# ggsave("../figures/p1_time_series.pdf", p_time_series,
#        width=5, height=3, units="in")

## ---- p11_stats
# Bootstrap samples parameter estimates
beta_hat_la = t(apply(
  x_hat_la, 1, function(x){ARp.beta.est(x, 2)$LA}))
beta_hat_ls = t(apply(
  x_hat_ls, 1, function(x){ARp.beta.est(x, 2)}$LS))

# Estimates of bias and variance of the parameters
beta_hat_la_mean = colMeans(beta_hat_la)
beta_hat_ls_mean = colMeans(beta_hat_ls)
beta_hat_la_bias = beta_hat_la_mean - beta_la
beta_hat_ls_bias = beta_hat_ls_mean - beta_ls
beta_hat_la_var = apply(beta_hat_la, 2, var)
beta_hat_ls_var = apply(beta_hat_ls, 2, var)

## ---- p1_2
# Residuals from the original data with the bootstrap estimators
bootstrap_estimator_resids_la = t(apply(
  cbind(x_hat_la, beta_hat_la),
  1,
  function(x){ARp.resid(x[1:100], x[101:102])}
))
bootstrap_estimator_resids_ls = t(apply(
  cbind(x_hat_ls, beta_hat_ls),
  1,
  function(x){ARp.resid(x[1:100], x[101:102])}
))

# Predictions of x_101 using bootstrap samples of the resids above
AR2.bootstrap.pred = function(x0, beta, e){
  return(sum(beta * rev(x0)) + e)
}

x0 = data3A$x[99:100]
x101_preds_la = t(apply(
  cbind(beta_hat_la, bootstrap_estimator_resids_la),
  1,
  function(x){AR2.bootstrap.pred(x0, x[1:2], x[3:100])}
))
x101_preds_ls = t(apply(
  cbind(beta_hat_ls, bootstrap_estimator_resids_ls),
  1,
  function(x){AR2.bootstrap.pred(x0, x[1:2], x[3:100])}
))

x101_mean_la = mean(x101_preds_la)
x101_mean_ls = mean(x101_preds_ls)
x101_pi_la = quantile(x101_preds_la, c(0.025, 0.975))
x101_pi_ls = quantile(x101_preds_ls, c(0.025, 0.975))

# Plot some of the predictions
nums2 = c(10, 50, 70)
p_preds = p_x
for(i in 1:length(nums2)){
  p_preds = p_preds +
    geom_point(
      data=tibble(x=x101_preds_la[nums[i], nums2[i]], t=101),
      aes(t, x, col="LA"), shape="x", size=2) +
    geom_point(
      data=tibble(x=x101_preds_ls[nums[i], nums2[i]], t=101),
      aes(t, x, col="LS"), shape="x", size=2) +
    labs(col="")
}
# ggsave("../figures/p1_preds.pdf", p_preds,
#        width=5, height=3, units="in")

## ---- save
save(beta_la, beta_ls,
     beta_hat_la_bias, beta_hat_ls_bias,
     beta_hat_la_var, beta_hat_ls_var,
     file="../data/p1.Rdata")
