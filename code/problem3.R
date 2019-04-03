setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- p3_prelim
# Load libraries and data
library(ggplot2)
library(tibble)
set.seed(123)
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


## ---- bs3c
library(gridExtra)
set.seed(56)
n <- length(data$u)
lambda_init = c(1,1)
B <- 1500
theta.hat.boot <- matrix(NA,nrow = B, ncol = 2)
for (i in seq(1,B)){
  index = sample(1:n,n,replace = T)
  theta.hat.boot[i,] = em(data$z[index], data$u[index], lambda_init, 1e-2)$lambda
}

lambda.df <- data.frame(lambda0 = theta.hat.boot[,1],lambda1 = theta.hat.boot[,2])
load(file="../data/p3.Rdata")

sd.lambda0 <- sd(lambda.df$lambda0)
sd.lambda1 <- sd(lambda.df$lambda1)


bias.lambda0 <- mean(lambda.df$lambda0) - lambda[1]
bias.lambda1 <- mean(lambda.df$lambda1) - lambda[2]

cor.lambda <- cor(lambda.df$lambda0, lambda.df$lambda1)

# Calculating the bias correction 
lambda0c = lambda[1]- bias.lambda0
lambda1c = lambda[2] - bias.lambda1 


hist.lambda0 <- ggplot(lambda.df) + 
  geom_histogram(aes(x = lambda0, y = ..density..),bins = 30)+ 
  geom_vline(xintercept = lambda[1], color = "firebrick",size=1) 

hist.lambda1 <- ggplot(lambda.df) + 
  geom_histogram(aes(x = lambda1, y = ..density..),bins = 30)+ 
  geom_vline(xintercept = lambda[2], color = "firebrick",size=1)

hist.lambda <- grid.arrange(hist.lambda0,hist.lambda1)

cat("Standard diviation \\lambda_0:", sd.lambda0)
cat("Standard diviation \\lambda_1:", sd.lambda1)
cat("Bias of \\hat{\\lambda}_0:", bias.lambda0)
cat("Bias of \\hat{\\lambda}_1:", bias.lambda1)
cat("Correlation between \\hat{\\lambda}_0 and \\hat{\\lambda}_1:",cor.lambda)
cat("Bias corrected \\lambda_0:", lambda0c)
cat("Bias corrected \\lambda_1:", lambda1c)


## ---- break
ggsave("../figures/hist_lambda.pdf", hist.lambda,
       width=5.5, height=5.5, units="in")
c3print <- c(sd.lambda0,sd.lambda1,bias.lambda0,bias.lambda1,cor.lambda,lambda0c,lambda1c)
save(file = "../data/variables/c3print.Rdata", c3print)


## ---- print3c
load(file = "../data/variables/c3print.Rdata")
cat("Standard diviation \\lambda_0:", c3print[1])
cat("Standard diviation \\lambda_1:", c3print[2])
cat("Bias of \\hat{\\lambda}_0:", c3print[3])
cat("Bias of \\hat{\\lambda}_1:", c3print[4])
cat("Correlation between \\hat{\\lambda}_0 and \\hat{\\lambda}_1:",c3print[5])
cat("Bias corrected \\lambda_0:", c3print[6])
cat("Bias corrected \\lambda_1:", c3print[7])
## ---- bs3d