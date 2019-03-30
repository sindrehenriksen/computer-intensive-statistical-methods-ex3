setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- p3_prelim
# Load libraries and data
library(ggplot2)
library(tibble)
data <- data.frame(
  z = read.table("../data/z.txt",header=F, col.names = "z"),
  u = read.table("../data/u.txt",header=F, col.names = "u")
)

## ---- p3_2
# M-step function (returns argmax of complete likelihood)
argmax_Q = function(z, u, lambda) {
  n = length(z)
  stopifnot(length(u) == n)
  stopifnot(length(lambda) == 2)

  lambda_0 = n * sum(
    u * z +
    (1 - u) * (1 / lambda[1] - z / (exp(lambda[1] * z) - 1))
  )^(-1)
  lambda_1 = n * sum(
    (1 - u) * z +
    u * (1 / lambda[2] - z / (exp(lambda[2] * z) - 1))
  )^(-1)

  return(c(lambda_0, lambda_1))
}

# Expecation maximation algorithm
em = function(z, u, lambda_initial, tol) {
  step = Inf
  lambda = lambda_initial
  i = 0
  steps = numeric()
  lambda_steps = list(lambda_initial)
  while(step > tol) {
    i = i + 1
    temp = lambda
    lambda = argmax_Q(z, u, lambda)
    step = sqrt(sum((lambda - temp)^2))
    steps[i] = step
    lambda_steps[[i + 1]] = lambda
  }
  print(paste0("Converged after ", i, " iterations"))
  return(
    list(lambda=lambda, lambda_steps=lambda_steps, steps=steps))
}

# Estimate lambda
lambda_init = c(1, 1)
res = em(data$z, data$u, lambda_init, 1e-2)
lambda = res$lambda

# Plot convergence
steps = res$steps
steps_df = tibble(iteration=(1:length(steps)), step=steps)
p_convergence = ggplot(steps_df, aes(x=iteration, y=step)) +
  geom_point() +
  labs(y="l2-norm of step size")
# ggsave("../figures/p3_convergence.pdf", p_convergence,
       # width=4, height=3, units="in")

lambda_steps = do.call("rbind", res$lambda_steps)
lambda_steps_df = tibble(iteration=0:length(steps),
                      lambda0=lambda_steps[,1],
                      lambda1=lambda_steps[,2])
p_lambda = ggplot(lambda_steps_df) +
  geom_point(aes(x=iteration, y=lambda0, col="lambda0"),
             alpha=0.7) +
  geom_point(aes(x=iteration, y=lambda1, col="lambda1"),
             alpha=0.7) +
  labs(y="Lambda values", col="")
# ggsave("../figures/p3_lambda.pdf", p_lambda,
       # width=5, height=3, units="in")

# ggplot() +
  # geom_histogram(data=tibble(y=z), aes(y, ..density..)) +
  # geom_histogram(data=tibble(y=data$z), aes(y, ..density..), fill="2")

## ---- p3_save
save(lambda, file="../data/p3.Rdata")
