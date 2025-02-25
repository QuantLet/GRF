# Remove all objects in R
rm(list = ls())

# Getting the path of your current open file ----
current_path = rstudioapi::getActiveDocumentContext()$path 

setwd(dirname(current_path ))

# Loading source files ----
source("support/packages_and_libraries.R") 
source("support/functions.R") 

ptm <- proc.time()
# Data initializing  ------------------------------------------------------
tau = c(0.5) #confidence level
sig = 1 #sigma in error
#widths = c(0.001,0.01, 0.04, 0.1, 0.2) # used for the plot along with for loop
width = 0.2 #eeta
c = 5 #constant
n1 = 1000 #data points for x1 
b = 500 #number of bootstraps for MBS
reps = 100 #repetitions of whole simulation
grids_x1 = 50 #grid points for CBs
grids_x1_list = c(15,20,30,50, 100)
set.seed(100)
node_size = c(5) #h
node_sizes = c(2,5,7,10,15) #bias controlling param of RFS
x2_fixed = c(0.3,0.5) #fix x2
alpha_sig = 0.05

rfs = list()
T_stats = list()
CIs = list()
CIs_std = list()
coverages = list()
coverages_std = list()
coverages_widths_PW = list()
coverages_widths_uniform = list()
coverages_PW_std = list()
coverages_uniform_std = list()
T_stat = list()
size_bandwidth_uni = list()
size_bandwidth_std_uni = list()

### here define calculation type 'size', 'power' or 'normal'
cal_type = 'normal'
## if power, provide a constant
theta_0s = seq(0.00001, 1, length.out = 5)
theta_0 = 0


# Data generation ----
#for (grids_x1 in grids_x1_list){
#for (node_size in node_sizes){
#for (theta_0 in theta_0s){
set.seed(123)
## Sample
X = get_x(n1, seed=42) #no replications for X because X is deterministic
theta_fun = function(X) theta_triangle(X, cal_type = cal_type, c = theta_0 )

#theta_true = theta_triangle(X, cal_type = cal_type, c = theta_0 ) #+ qnorm(tau)*sig ##this is used for the case of quantile reg
#Y = get_y(X, zero, sig, NULL,reps) #for size
Y = get_y(X, theta_fun, sig, NULL,reps)
n = nrow(X)

## Test set
X_test = get_x(grids_x1, seed = 123)
theta_true_test = theta_triangle(X_test,  cal_type = cal_type, c = theta_0)# + qnorm(tau)*sig
grids = nrow(X_test)

## Fitting the forest ----
rf= fit_forest(X,Y,tau = tau,node_size = node_size)
preds = predict_forest(rf,X_test)
theta_hat_test = lapply(1:reps, function(j) preds[[j]]$predictions)
w_test = weights_forest(rf,X_test)
mat <- do.call("cbind",theta_hat_test)
theta_hat_expected = rowMeans(mat)


## Calculations for quantile forest ----
# kde = lapply(1:reps, function (j) density(Y[,j], n=nrow(X))) #estimation of the density
# f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat_test[[1]]))[2]))
# V_hat = 1/f_Y
# psi = lapply(1:reps, function(k) t(sapply(1:nrow(X_test), function(j) tau - (Y[,k]<=theta_hat_test[[k]][j]))))
# H_hat = sapply(1:reps, function(j) n* apply(w_test[[j]]*psi[[j]],1,var))
# sigma_hat = (f_Y^(-2)*H_hat)^(1/2)
# set.seed(5)
# e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))
# T_stat = lapply(1:reps, function(k) sapply(1:b, function(j)
#   (-(H_hat[,k]^(-1/2)*w_test[[k]]*psi[[k]])%*% e_multipliers[[k]][[j]])@x))
# #T_stat = test_stat(X,Y,X_test,Y_test,theta_hat_test,w_test, tau= tau, b)
# T_stat_abs = lapply(1:reps, function(j) abs((T_stat[[j]])))


## Calculation for regression forest
V_hat = -1
psi = lapply(1:reps, function(k) t(sapply(1:nrow(X_test), function(j) (Y[,k]-theta_hat_test[[k]][j]))))
H_hat = sapply(1:reps, function(j) n* apply(w_test[[j]]*psi[[j]],1,var))
#H_hat = sapply(1:reps, function(j) w_test[[j]]*psi[[j]])

#sigma_hat = (H_hat)^(1/2)

## alternative formulation for forest variance
var_hat = sapply(1:reps, function(j) preds[[j]]$variance.estimates)
sigma_hat = (var_hat)^(1/2)


e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))
T_stat = lapply(1:reps, function(k) sapply(1:b, function(j)
  ((sigma_hat[,k]^(-1)*w_test[[k]]*psi[[k]])%*% e_multipliers[[k]][[j]])))
T_stat_abs = lapply(1:reps, function(j) abs((T_stat[[j]])))



## just for simplicity, renaming test sets as original sets
# X = X_test
# Y = Y_test
# theta_hat = theta_hat_test
# theta_true = theta_true_test
# w = w_test

## Confidence interval with test stat ----
### alpha quantile of test statistic 
q_star = lapply(1:reps, function(j) rowQuantiles(T_stat_abs[[j]] , probs= c(1-alpha_sig)))
CI = confidence_interval(theta_hat_test, q_star, sigma_hat,reps)
## alpha quantile for standard normal
set.seed(123)
q_norm = lapply(1:reps, function(j) qnorm(1-alpha_sig/2))
CI_std = confidence_interval(theta_hat_test, q_norm, sigma_hat,reps)


## Calculating the coverage ----
### Coverage of true theta
coverage_bootstrap_pw = coverage(theta_true_test, CI, grids, reps)
coverage_std_pw = coverage(theta_true_test, CI_std, grids, reps)
### Coverage of expected theta
coverage_expected_pw = coverage(theta_hat_expected, CI, grids, reps)
#coverage_expected_std = coverage(theta_hat_expected, CI_std, grids, reps)
#T_stats[[as.character(n)]]  = T_stat
#CIs[[as.character(n)]]  = CI
#coverages[[as.character(n)]] = coverage_
#coverages_std[[as.character(n)]] = coverage_std


## Calculation of uniform confidence bands ----
uniform_T_stat = lapply(1:reps, function(j) T_stat_abs[[j]]) 
uniform_T_max = lapply(1:reps, function(j) apply(uniform_T_stat[[j]], 2, max)) # max of t_stats 
uniform_q_star = lapply(1:reps, function(j) quantile(uniform_T_max[[j]], 1-alpha_sig)) # quantile of max_t_stat

uniform_CI = confidence_interval(theta_hat_test, uniform_q_star, sigma_hat, reps)


coverage_bootstrap_uniform = coverage_uniform(theta_true_test, uniform_CI, grids,reps)
coverage_expected_uniform = coverage_uniform(theta_hat_expected, uniform_CI, grids,reps)


## calculation of asymptotic uniform coverage
std_T_stat = lapply(1:reps,function(k)  sapply(1:b, function(j) abs(rnorm(grids))))
std_T_max = lapply(1:reps, function(j) apply(std_T_stat[[j]], 2, max))
std_q_star = lapply(1:reps, function(j) quantile(std_T_max[[j]], 1-alpha_sig)) # quantile of max_t_stat
std_CI = confidence_interval(theta_hat_test, std_q_star, sigma_hat, reps)
coverage_std_uniform = coverage_uniform(theta_true_test, std_CI, grids,reps)


#coverages_widths_PW[[as.character(grids_x1)]] = coverage_bootstrap_pw
#coverages_widths_uniform[[as.character(grids_x1)]] = coverage_bootstrap_uniform

#coverages_widths_PW[[as.character(node_size)]] = coverage_bootstrap_pw
#coverages_widths_uniform[[as.character(node_size)]] = coverage_bootstrap_uniform

#coverages_widths_PW[[as.character(theta_0)]] = coverage_bootstrap_pw
#coverages_widths_uniform[[as.character(theta_0)]] = coverage_bootstrap_uniform
#coverages_PW_std[[as.character(theta_0)]] = coverage_std_pw
#coverages_uniform_std[[as.character(theta_0)]] = coverage_std_uniform

## calculating size
# Size is defined as prob of type 1 error (false positive) given true theta is 0
# 
size_pw = 1-(coverage_bootstrap_pw/100)
size_std_pw = 1-(coverage_std_pw/100)
size_uni = 1-(coverage_bootstrap_uniform/100)
size_std_uni = 1-(coverage_std_uniform/100)

size_bandwidth_uni[[as.character(node_size)]] = size_uni
size_bandwidth_std_uni[[as.character(node_size)]] = size_std_uni

#}

## MSE for a fixed test point 
mse = rowMeans((do.call("cbind",theta_hat_test)-theta_true_test)^2)[1]


print(proc.time() - ptm)

# removing undesirable variables ----
rm(ptm)

# Bias-Variance Analysis
df <- as.data.frame(do.call(cbind, theta_hat_test))
empirical <- apply(df, 1, sd, na.rm=TRUE)
analytical <- apply(sigma_hat,1, mean)
var_ratio <- mean(empirical/analytical)
bias <- mean(apply(sapply(1:reps, function(j) abs(theta_true_test - theta_hat_test[[j]])),1,mean))


# Print size ----

# Print variables in one line
cat(    "n =", n1, 
    "grid =", grids_x1, 
    "h =", node_size, 
    "sig =", sig, "\n")

cat("Size (Bootstrap Pointwise):", size_pw, "\n")
cat("Size (Standard Pointwise):", size_std_pw, "\n")
cat("Size (Bootstrap Uniform):", size_uni, "\n")
cat("Size (Standard Uniform):", size_std_uni, "\n")


# Print coverages ----
cat("Coverage (Bootstrap Pointwise):", coverage_bootstrap_pw, "\n")
cat("Coverage (Standard Pointwise):", coverage_std_pw, "\n")
cat("Coverage (Bootstrap Uniform):", coverage_bootstrap_uniform, "\n")
cat("Coverage (Standard Uniform):", coverage_std_uniform, "\n")





