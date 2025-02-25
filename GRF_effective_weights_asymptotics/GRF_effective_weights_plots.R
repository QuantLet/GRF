#install.packages("matrixStats") ----
rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
library(matrixStats)
set.seed(42)


# Getting the path of your current open file ----
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Function definition for data --------------------------------------------
get_x = function(n,seed){
  set.seed(seed)
  X= expand.grid(X1 = seq(-0.5, 0.5, length.out = n))
  #X= expand.grid(X1 = seq(-0.5, 0.5, length.out = n), X2 = seq(0, 1, 0.1))
  return(X)
}

theta_triangle = function(x, width){
  #pmax(1 - abs((x[[1]]) / width), 0)
  x[[1]]+4*x[[1]]^3
} 

get_y = function(X, theta, sigma, seed=NULL, reps = 1){
  set.seed(NULL)
  n = nrow(X)
  return(replicate(reps,(theta(X) + rnorm(n, 0, sigma))))
} 

# Data initializing  ------------------------------------------------------
tau = c(0.5)
sig = 0.1
width = 0.2
reps = 10 #repititions
grids = 100 #grid points
n_vals = c(500, 2000, 3000)  # Sample sizes

set.seed(100)

node_size = 3
errors <- c()
for (n in n_vals) {
  
  ## Estimation ---------------
  set.seed(100)
  X = get_x(n, seed=42) #no replications for X because X is deterministic
  theta_fun = function(X) theta_triangle(X, width)
  theta_true = theta_triangle(X, width) + qnorm(tau)*sig
  Y = get_y(X, theta_fun, sig, NULL,reps)
  rand_for =  function(j)  grf::quantile_forest( X ,data.matrix(Y[,j]), 
                                                 quantiles = tau, min.node.size = node_size)
  rf = lapply(1:reps, rand_for)
  w = sapply(1:reps,function(j) get_forest_weights(rf[[j]]))
  objective_fun = function(theta,Y,alpha) 
    sum(((Y-theta)) * as.matrix(alpha) * (tau -   (Y <= theta)))
  theta_hat = lapply(1:reps, function(j) sapply(1:nrow(X), function(k)
    optimize(f=objective_fun,interval = c(0,1),
             tol = 0.0001, Y=Y[,j], alpha=w[[j]][k,])[1]))
  
  
  ## Calculations for test set ----
  X_test = get_x(grids, seed=50)
  Y_test = get_y(X_test, theta_fun, sig, NULL,reps)
  theta_hat_test = lapply(1:reps, function(j) approx(X$X1,unlist(theta_hat[[j]]),
                                                     xout = X_test$X)[["y"]])
  rand_for =  function(j)  grf::quantile_forest( X_test ,data.matrix(Y_test[,j]),
                                                 quantile = tau, min.node.size = node_size)
  rf = lapply(1:reps, rand_for)
  #w = sapply(1:reps,function(j) get_forest_weights(rf[[j]]))
  w_test = lapply(1:reps,function(j) get_forest_weights(rf[[j]]))
  theta_true_test = theta_triangle(X_test, width) + qnorm(tau)*sig
  
  ## just for simplicity, renaming test sets as original sets
  
  
  kde = lapply(1:reps, function (j) density(Y_test[,j], n=n)) #estimation of the density
  f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat_test[[1]]))[2]))
  
  epsilon_tilde <- sapply(1:reps, function(j) -f_Y^(-1)*(tau - (Y_test[,j] <= unlist(theta_hat_test[[j]]))))
  
  theta_tilde_test <- lapply(1:reps, function(k) theta_true_test +
                               (w_test[[k]] %*% ((f_Y[,k]^(-1)) * (tau - (Y_test[,k] <= unlist(theta_hat_test[[k]]))))))
  
  sup_diff <- sapply(1:reps, function(j) max(abs(theta_hat_test[[j]] - theta_tilde_test[[j]])))
  average_error <- mean(sup_diff)
  
  errors <- c(errors, average_error)
  
  # Plot for one rep
  
  rep <-1 #fixed replication
  
  plot(X_test[, "X1"], Y_test[,rep], pch = 20, cex = 0.5, col = "gray", main = glue("n = {n}, Realization {rep}"),
       xlab = 'X', ylab = 'Y', ylim = c(-2.5,2.5))
  lines(X_test[, "X1"], theta_true_test, col = "blue", lty = 2)
  lines(X_test[, "X1"], theta_tilde_test[[rep]], col = "red")
  lines(X_test[, "X1"], theta_hat_test[[rep]], col = "green")
  legend("topright", legend = c("True Theta", "Theta Tilde", "Theta Hat"), 
         col = c("blue", "red", "green"), lty = c(2, 1, 1), cex = 0.8)
  
}

print(errors)

