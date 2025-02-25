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
  pmax(1 - abs((x[[1]]) / width), 0)
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
c = 1
b = 100 #number of bootstraps for MBS
reps = 10 #repititions for confidence interval
grids = 100 #grid points for CBs
n= 500
set.seed(100)

T_stats = list()
CIs = list()
CIs_std = list()
coverages = list()
coverages_std = list()
node_size = 3


## Estimation ---------------
T_stat = list()
set.seed(100)
ptm <- proc.time()
X = get_x(n, seed=42) #no replications for X because X is deterministic
theta_fun = function(X) theta_triangle(X, width)
theta_true = theta_triangle(X, width) + qnorm(tau)*sig
Y = get_y(X, theta_fun, sig, NULL,reps)
rand_for =  function(j)  grf::quantile_forest( X ,data.matrix(Y[,j]), 
                              quantiles = tau, min.node.size = node_size)
rf = lapply(1:reps, rand_for)
w = sapply(1:reps,function(j) get_sample_weights(rf[[j]]))
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
w = sapply(1:reps,function(j) get_sample_weights(rf[[j]]))
w_test = lapply(1:reps,function(j) get_sample_weights(rf[[j]]))
theta_true_test = theta_triangle(X_test, width) + qnorm(tau)*sig

## just for simplicity, renaming test sets as original sets
X = X_test
Y = Y_test
theta_hat = theta_hat_test
theta_true = theta_true_test
w = w_test

kde = lapply(1:reps, function (j) density(Y[,j], n=n)) #estimation of the density
f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat[[1]]))[2]))
V_hat = 1/f_Y
H_hat = sapply(1:reps, function(j) sapply(1:nrow(X), 
                                          function(k)  n *((var(w[[j]][k,]*(tau - (Y[,j] <= unlist(theta_hat[[j]]))))))))
sigma_hat = (f_Y^(-2)*H_hat)^(1/2)
e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(nrow(X), 0, 1)))

T_stat = lapply(1:reps, function(k) sapply(1:b, function(j) 
  (w[[k]] %*% ((H_hat[,k]^(-1/2)) * (tau - (Y[,k] <= unlist(theta_hat[[k]])))  * e_multipliers[[k]][[j]]))@x))


## Confidence interval with test stat
alpha_sig = 0.05
T_stat_abs = lapply(1:reps, function(j) abs(t(T_stat[[j]])))
q_star = lapply(1:reps, function(j) colQuantiles(T_stat_abs[[j]] , probs= c(1-alpha_sig)))
CI = lapply(1:reps, function(j) list(unlist(theta_hat[[j]])-(q_star[[j]]*sigma_hat[j]),
                                     unlist(theta_hat[[j]])+(q_star[[j]]*sigma_hat[j])))
print(proc.time() - ptm)

## Calculating the coverage  ----

T_stats[[as.character(n)]]  = T_stat
CIs[[as.character(n)]]  = CI

q_norm = qnorm(1-alpha_sig)
CI_std = lapply(1:reps, function(j) list(unlist(theta_hat[[j]])-(q_norm*sigma_hat[j]),
                                         unlist(theta_hat[[j]])+(q_norm*sigma_hat[j])))
coverage = mean(sapply(1:reps, function(k)
  sum((theta_true > CI[[k]][[1]]) & (theta_true < CI[[k]][[2]]))/nrow(X)))*100
coverage_std = mean(sapply(1:reps, function(k)
  sum(theta_true > CI_std[[k]][[1]] & theta_true < CI_std[[k]][[2]])/nrow(X)))*100

T_stats[[as.character(n)]]  = T_stat
CIs[[as.character(n)]]  = CI
coverages[[as.character(n)]] = coverage
coverages_std[[as.character(n)]] = coverage_std


## Uniform confidence bands ----

#grid_T_stat = T_stat_abs[[1]][,seq(1, nrow(X), length.out = grids) ] #t stats for CB
grid_T_stat = T_stat_abs[[1]] #t stats for CB
grid_T_max = apply(grid_T_stat, 1, max)    
grid_q_star = quantile(grid_T_max, 1-alpha_sig)
pd = data.frame(X1=X$X1, sigma = sigma_hat[[1]],theta_hat =  unlist(theta_hat[[1]]),
                theta_true = theta_true, CI_L = CI[[1]][[1]], CI_U = CI[[1]][[2]])
pd$grid_CI_L = pd$theta_hat-(grid_q_star*pd$sigma)
pd$grid_CI_U = pd$theta_hat+(grid_q_star*pd$sigma)


# Plotting ----
png(file = glue('CI_bands_','n{formatC(as.integer(n*10), width=4, flag="0")}_',
                'tau{formatC(tau*10, width=3 ,flag="0")}_',
                'sig{formatC(sig*10, width=3 ,flag="0")}_',
                'grids{formatC(as.integer(grids) ,flag="0")}_',
                '.png'),
    width=1500, height=1500)

plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), 
     ylim=range(c(pd$theta_true,pd$theta_hat, pd$CI_L,pd$CI_U, pd$grid_CI_U,pd$grid_CI_L)))
points(pd$X1, pd$theta_true,
       col='red', main="Confidence intervals", pch=19, cex=2)
points(pd$X1, pd$theta_hat, col='blue', pch=19, cex=2)
lines(pd$X1, pd$CI_L,      col='black',pch = 19,type = "b", lty = 2, cex=0.8)
lines(pd$X1, pd$CI_U,      col='black',pch = 19,type = "b", lty = 2, cex=0.8)
lines(pd$X1, pd$grid_CI_L, col='magenta',pch = 19, lty = 2, cex=2)
lines(pd$X1, pd$grid_CI_U, col='magenta',pch = 19, lty = 2, cex=2)
dev.off()

# # plot with average confidence intervals (cancelled) ----
# avg_ci_lower = 0
# avg_ci_upper = 0
# avg_theta_hat = 0
# for (i in seq(1, reps,1)){
#   avg_ci_lower = avg_ci_lower + CI[[i]][[1]] 
#   avg_ci_upper = avg_ci_upper + CI[[i]][[2]]
#   avg_theta_hat = avg_theta_hat + unlist(theta_hat[[i]])
# }
# avg_ci_lower = avg_ci_lower/reps
# avg_ci_upper = avg_ci_upper/reps
# avg_theta_hat = avg_theta_hat/reps
# 
#   png(file = glue('CI_averaged_','n{formatC(as.integer(n), width=4, flag="0")}_',
#                   'sig{formatC(sig*10, width=3 ,flag="0")}_',
#                   'reps{formatC(as.integer(reps) ,flag="0")}_',
#                   '.png'),
#       width=1500, height=1500)
#   
#   pd = data.frame(X1=X$X1, theta_hat =  avg_theta_hat,
#                   theta_true = theta_true, CI_L = avg_ci_lower, CI_U = avg_ci_upper)
#   pd = pd[order(X$X1),] #ordering X wrt X1
#   plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), 
#        ylim= range(pd$theta_true, pd$theta_hat, pd$CI_L, pd$CI_U))
#   points(pd$X1, pd$theta_true ,
#          ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
#          col='red', main="Confidence intervals", pch=19, cex=2)
#   points(pd$X1, pd$theta_hat, col='blue', pch=19, cex=2)
#   lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
#   lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
#   dev.off()
# 
# ## plot with multiple confidence intervals 
# for (num in c(1,1)){
#     png(file = glue('CI_','n{formatC(as.integer(n), width=4, flag="0")}_',
#                     'sig{formatC(sig*10, width=3 ,flag="0")}_',
#                     'reps{formatC(as.integer(num) ,flag="0")}_',
#                     '.png'),
#         width=1500, height=1500)
#     
#     plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5),
#          ylim= range(pd$theta_true, pd$theta_hat, pd$CI_L, pd$CI_U))
#     for (i in 1:num){
#       pd = data.frame(X1=X$X1, theta_hat =  unlist(theta_hat[[i]]),
#                       theta_true = theta_true, CI_L = CI[[i]][[1]], CI_U = CI[[i]][[2]] )
#       pd = pd[order(X$X1),] #ordering X wrt X1
#       points(pd$X1, pd$theta_true ,
#              ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
#              col='red', main="Confidence intervals", pch=19, cex=2)
#       points(pd$X1, pd$theta_hat, col='blue', pch=19, cex=2)
#       lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
#       lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
#     }
#     dev.off()
#   }