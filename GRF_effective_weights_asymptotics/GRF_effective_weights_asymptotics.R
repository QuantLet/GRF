rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
#install.packages('ks')
library(ks)
set.seed(42)


get_x = function(n, c){
  X = matrix(sort(runif(n,-c, c)), nrow = n)
  colnames(X) = 'X1'
  return(X)
}
get_x_grid = function(n, c){
  X = matrix(seq(-c, c, length.out = n), nrow = n)
  colnames(X) = 'X1'
  return(X)
}
polynomial = function(x, p){
  X = matrix(1, nrow = nrow(x), ncol = ncol(x) * p + 1)
  cols = c('const')
  for (j in 1:ncol(x)){
    for (i in 1:p) {
      X[, (j - 1) * ncol(x) + i + 1] = x[, j] ** i
      cols = c(cols, glue("{colnames(x)[j]}^{i}"))
    }
  }
  colnames(X) = cols
  return(X)
}
theta_polynomial = function(x, p, beta){
  return(polynomial(x, p) %*% beta)
}

get_y = function(X, theta, sigma, seed=NULL){
  set.seed(NULL)
  n = nrow(X)
  return(theta(X) + rnorm(n, 0, sigma))
}

tau = c(0.5)
sig = 1
width = 0.2
beta = c(0, 1, 0, 4)
p = 3
c = 1
grids = 50
reps = 100



rfs = list()
rmse_results = list()
#n=500
#rep = 1
for (n in c(500, 1000, 2000)){
  sup_diffs = numeric(reps)   # Store sup-norm results for current n
  set.seed(123)
  X_test = get_x(grids, c)  # New test set
  for (rep in 1:reps) {  
  X = get_x(n, c)
  theta = function(X) theta_polynomial(X, p, beta)
  Y = get_y(X, theta, sig, 42)
  rfs[[as.character(n)]] = rf = grf::quantile_forest(X, Y, quantiles = tau, seed = 42)
  
   theta_test = theta(X_test)
   theta_hat = predict(rf, X_test)$predictions
  alpha = get_forest_weights(rf, X_test)  # Compute weights for test set
  
  kde_fit<- ks::kde(Y)  # Kernel density estimation
  f_Y_theta <- - (predict(kde_fit, x = theta_hat) ) ^ (-1)
  psi_matrix <- sapply(theta_hat, function(theta_t) tau - (Y <= theta_t))
  
  epsilon_tilde <- sapply(seq_along(theta_hat), function(i) -predict(kde_fit, x= theta_hat[i])^(-1) * psi_matrix[,i])
  
  theta_tilde = (theta_test + (alpha %*% epsilon_tilde)[1])
  
  #Compute sup-norm of the absolute difference
  sup_diffs[rep] = max(abs(theta_hat - theta_tilde))
  }
  
  # Compute RMSE comparing sup-norms to true theta
  true_theta_values = theta(X)
  rmse = mean(sup_diffs) 
  rmse_results[[as.character(n)]] = rmse
  
}

# Print RMSE results
print(rmse_results)

