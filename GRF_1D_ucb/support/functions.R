#################################
# This file contains all the functions required for running the
# simulations for constructing uniform confidence bands for GRF estimates


# Functions definition for data --------------------------------------------
get_x = function(n,seed=NULL){
  set.seed(seed)
  X= expand.grid(X1 = seq(-0.5, 0.5, length.out = n))#, X2 = seq(0.1, 1, 0.1))
  return(X)
}

theta_triangle = function(x, cal_type, c = NULL){
  
  if(cal_type == 'size'){
    theta = 0
  } 
  
  else {
    
    ## uncomment your desired function
    #theta = abs(x[[1]])+0.5
    #n = length(x[[1]])
    #theta = seq(0, 1, length.out = n)
    theta = pmax(1 - abs((x[[1]]) / 0.2), 0)
    #theta =  sin(x[[1]]*3*pi)
    #theta =  sin(x[[1]]*8*pi)
    #theta =  sin(x[[1]]*5*pi)
    #theta =  sin(x[[1]]*3)
    #theta =  sin(x[[1]]*5)
    #theta =  sin(x[[1]]*8)
    #theta = sin(x[[1]])
    #theta = 1+x[[1]]+2*x[[1]]^2+3*x[[1]]^3
  }
  return(theta)
} 

zero = function(x){
  return(0)
}

get_y = function(X, theta, sigma, seed=10, reps = 1){
  set.seed(10)
  n = nrow(X)
  return(replicate(reps,(theta(X) + rnorm(n, 0, sigma))))
} 

# Functions for effective weights and RF estimated parameter for test points ----

fit_forest = function(X,Y,tau = 0.5,node_size = 5) { 
  #rand_for =  function(j)  grf::quantile_forest( X ,data.matrix(Y[,j]), quantile = tau, min.node.size = node_size)
  rand_for =  function(j)  grf::regression_forest( X ,data.matrix(Y[,j]), min.node.size = node_size, num.trees = 1000)
  reps = dim(Y)[2]
  rf = lapply(1:reps, rand_for)
  return(rf)
}

predict_forest = function(rf, X_test){
  reps = length(rf)
  preds= lapply(1:reps, function(j) predict(rf[[j]],X_test, estimate.variance = TRUE))
  return(preds)
}

weights_forest = function(rf, X_test){
  reps = length(rf)
  w_test = lapply(1:reps,function(j) as.matrix(get_forest_weights(rf[[j]], newdata = X_test)))
  return(w_test)
}

# Functions for calculation of test statistic ----
test_stat = function(X,Y,X_test,Y_test,theta_hat_test,w_test, tau= 0.5, b){
reps = dim(Y)[2]
n = nrow(X)
psi = lapply(1:reps, function(k) t(sapply(1:nrow(theta_hat_test[[k]]), function(j) tau - (Y[,k]<=theta_hat_test[[k]][j]))))
kde = lapply(1:reps, function (j) density(Y_test[,j], n=n)) #estimation of the density
f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat_test[[1]]))[2]))
V_hat = 1/f_Y
H_hat = sapply(1:reps, function(j) n* apply(w_test[[j]]*psi[[j]],1,var))
sigma_hat = (f_Y^(-2)*H_hat)^(1/2)
set.seed(5)
e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))
T_stat = lapply(1:reps, function(k) sapply(1:b, function(j)
  (-(H_hat[,k]^(-1/2)*w_test[[k]]*psi[[k]])%*% e_multipliers[[k]][[j]])@x))
return(T_stat)
}
  
  
  
# Functions for CI ----
confidence_interval = function(theta, q, sigma,reps){
CI = lapply(1:reps, function(j) list(unlist(theta[[j]])-(q[[j]]*sigma[,j]),
                                     unlist(theta[[j]])+(q[[j]]*sigma[,j])))
return(CI)
}
# Functions for coverage ----
coverage = function (theta, CI, grids, reps){
  coverage = mean(sapply(1:reps, function(k)
    sum((theta > CI[[k]][[1]]) & (theta < CI[[k]][[2]]))/grids))*100
  return(coverage)
}

# Function for size of uniform inference ----
coverage_uniform = function(theta, CI, grids, reps){
  coverage_uni = mean(sapply(1:reps, function(k) sum((theta > CI[[k]][[1]]) & (theta < CI[[k]][[2]]))==grids))*100
  return(coverage_uni)
  }