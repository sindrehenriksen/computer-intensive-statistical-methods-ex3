setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- p1_prelim
# Load libraries and data
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
p_x = ggplot(data) + geom_line(aes(t, x, col="Data"))
nums = c(300, 600, 1000)
p_time_series = p_x
for(i in 1:length(nums)){
  p_time_series = p_time_series +
    geom_line(
      data=tibble(x=x_hat_la[nums[i],], t=1:100),
      (aes(t, x, col="LA"))) +
    geom_line(
      data=tibble(x=x_hat_ls[nums[i],], t=1:100),
      (aes(t, x, col="LS"))) +
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
# Predictions of x_101 using random resid for all betas
AR2.bootstrap.pred = function(x0, beta, e){
  return(sum(beta * rev(x0)) + sample(e, 1))
}

x0 = data3A$x[99:100]
x101_preds_la = apply(
  beta_hat_la,
  1,
  function(x){AR2.bootstrap.pred(x0, x, epsilon_la)}
)
x101_preds_ls = apply(
  beta_hat_ls,
  1,
  function(x){AR2.bootstrap.pred(x0, x, epsilon_ls)}
)

x101_mean_la = mean(x101_preds_la)
x101_mean_ls = mean(x101_preds_ls)
x101_pi_la = quantile(x101_preds_la, c(0.025, 0.975))
x101_pi_ls = quantile(x101_preds_ls, c(0.025, 0.975))

# Plot 95% prediction intervals
p_preds = p_x
p_preds = p_preds +
  geom_ribbon(
    data=tibble(xmin=c(data3A$x[100], x101_pi_la[1]),
                xmax=c(data3A$x[100], x101_pi_la[2]),
                t=c(100, 101)),
    aes(x=t, ymin=xmin, ymax=xmax, fill="LA"), alpha=0.3) +
  geom_ribbon(
    data=tibble(xmin=c(data3A$x[100], x101_pi_ls[1]),
                xmax=c(data3A$x[100], x101_pi_ls[2]),
                t=c(100, 101)),
    aes(x=t, ymin=xmin, ymax=xmax, fill="LS"), alpha=0.3) +
  scale_fill_manual("", values=c(3, 4)) +
  labs(col="")
# ggsave("../figures/p1_preds.pdf", p_preds,
       # width=5, height=3, units="in")

## ---- p1_save
save(beta_la, beta_ls,
     beta_hat_la_bias, beta_hat_ls_bias,
     beta_hat_la_var, beta_hat_ls_var,
     x101_pi_la, x101_pi_ls,
     file="../data/p1.Rdata")
