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
  while(step > tol) {
    i = i + 1
    temp = lambda
    lambda = argmax_Q(z, u, lambda)
    step = sqrt(sum((lambda - temp)^2))
    steps[i] = step
  }
  print(paste0("Converged after ", i, " iterations"))
  return(list(lambda=lambda, steps=steps))
}

# Estimate lambda
lambda_init = c(1, 1)
res = em(data$z, data$u, lambda_init, 1e-2)
lambda = res$lambda

# Plot convergence
steps = res$steps
steps.df = tibble(iteration=(1:length(steps)), step=steps)
p_convergence = ggplot(steps.df, aes(x=iteration, y=step)) +
  geom_point() +
  labs(y="l2-norm of step size")
ggsave("../figures/p2_convergence.pdf", p_convergence,
       width=5, height=3, units="in")

# ggplot() +
  # geom_histogram(data=tibble(y=z), aes(y, ..density..)) +
  # geom_histogram(data=tibble(y=data$z), aes(y, ..density..), fill="2")
