# Packages installation ----
#install.packages("matrixStats")
#install.packages("matrixStats")
#install.packages('pracma')
#install.packages('interp')

# Packages loading ----
rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
library(matrixStats)
set.seed(42)

rm(list=ls())
library(grf)
library(ggplot2)
library(glue)
library(rstudioapi)
library(parallel)
library(matrixStats)
library(interp)
set.seed(42)

# Getting the path of your current open file ----
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Function definition for data --------------------------------------------
get_x = function(n,seed){
  set.seed(seed)
  X= expand.grid(X1 = seq(-0.5, 0.5, length.out = n), X2 = seq(0.1, 1, 0.1))
  return(X)
}

theta_triangle = function(x, width){
  pmax(1 - abs((x[[1]]) / width), 0)
} 

get_y = function(X, theta, sigma, seed=10, reps = 1){
  set.seed(10)
  n = nrow(X)
  return(replicate(reps,(theta(X) + rnorm(n, 0, sigma))))
} 

# Data initializing  ------------------------------------------------------
tau = c(0.5)
sig = 0.1
width = 0.2
c = 1
n1 = 200 #data points for x1 
b = 50 #number of bootstraps for MBS
reps = 20 #repititions of whole simulation
grids_x1 = 50 #grid points for CBs
set.seed(100)
node_size = 5 #bias controlling param of RFS
x2_fixed = c(0.3,0.5) #fix x2


rfs = list()
T_stats = list()
CIs = list()
CIs_std = list()
coverages = list()
coverages_std = list()

# Computation ----
T_stat = list()
set.seed(100)
ptm <- proc.time()
X = get_x(n1, seed=42) #no replications for X because X is deterministic
X = X[order(X$X1),] #ordering X wrt X1
n = nrow(X)
theta_fun = function(X) theta_triangle(X, width)
theta_true = theta_triangle(X, width) + qnorm(tau)*sig
Y = get_y(X, theta_fun, sig, NULL,reps)
rand_for =  function(j)  grf::quantile_forest( X ,data.matrix(Y[,j]), quantile = tau, min.node.size = node_size)
rf = lapply(1:reps, rand_for)
w = sapply(1:reps,function(j) get_forest_weights(rf[[j]]))


## Calculations for test set ----
X_test = expand.grid(X1 = seq(-0.5,0.5,length.out=grids_x1), X2 = x2_fixed)
Y_test = get_y(X_test, theta_fun, sig, NULL,reps)
grids = nrow(X_test)
theta_true_test = theta_triangle(X_test, width) + qnorm(tau)*sig
theta_hat_test = lapply(1:reps, function(j) predict(rf[[j]], newdata = X_test)$predictions)
w_test = lapply(1:reps,function(j) get_forest_weights(rf[[j]], newdata = X_test))

kde = lapply(1:reps, function (j) density(Y_test[,j], n=nrow(X))) #estimation of the density
f_Y= sapply(1:reps, function(j) unlist(approx(kde[[j]][["x"]], kde[[j]][["y"]], xout = c(theta_hat_test[[1]]))[2]))
V_hat = 1/f_Y
psi = lapply(1:reps, function(k) t(sapply(1:nrow(theta_hat_test[[k]]), function(j) tau - (Y[,k]<=theta_hat_test[[k]][j]))))
H_hat = sapply(1:reps, function(j) n* apply(w_test[[j]]*psi[[j]],1,var))
sigma_hat = (f_Y^(-2)*H_hat)^(1/2)
set.seed(5)
e_multipliers = lapply(1:reps, function(j) lapply(1:b, function(j) rnorm(n, 0, 1)))

T_stat = lapply(1:reps, function(k) sapply(1:b, function(j)
  (-(H_hat[,k]^(-1/2)*w_test[[k]]*psi[[k]])%*% e_multipliers[[k]][[j]])@x))

## just for simplicity, renaming test sets as original sets
X = X_test
Y = Y_test
theta_hat = theta_hat_test
theta_true = theta_true_test
w = w_test
## Confidence interval with test stat ----
alpha_sig = 0.05
T_stat_abs = lapply(1:reps, function(j) abs((T_stat[[j]])))
q_star = lapply(1:reps, function(j) colQuantiles(T_stat_abs[[j]] , probs= c(1-alpha_sig)))
CI = lapply(1:reps, function(j) list(unlist(theta_hat[[j]])-(q_star[[j]]*sigma_hat[j]),
                                     unlist(theta_hat[[j]])+(q_star[[j]]*sigma_hat[j])))

q_norm = qnorm(1-alpha_sig)
CI_std = lapply(1:reps, function(j) list(unlist(theta_hat[[j]])-(q_norm*sigma_hat[j]),
                                         unlist(theta_hat[[j]])+(q_norm*sigma_hat[j])))
print(proc.time() - ptm)


## Calculating the coverage ----
### Coverage of true theta
coverage_true = mean(sapply(1:reps, function(k)
  sum((theta_true > CI[[k]][[1]]) & (theta_true < CI[[k]][[2]]))/nrow(X)))*100
coverage_true_std = mean(sapply(1:reps, function(k)
  sum(theta_true > CI_std[[k]][[1]] & theta_true < CI_std[[k]][[2]])/nrow(X)))*100

### Coverage of expected theta
mat <- do.call("cbind",theta_hat)
theta_hat_expected = rowMeans(mat)
coverage_expected = mean(sapply(1:reps, function(k)
  sum((theta_hat_expected > CI[[k]][[1]]) & (theta_hat_expected < CI[[k]][[2]]))/nrow(X)))*100
coverage_expected_std = mean(sapply(1:reps, function(k)
  sum((theta_hat_expected > CI_std[[k]][[1]]) & (theta_hat_expected < CI_std[[k]][[2]]))/nrow(X)))*100

T_stats[[as.character(n)]]  = T_stat
CIs[[as.character(n)]]  = CI
coverages[[as.character(n)]] = coverage_true
coverages_std[[as.character(n)]] = coverage_true_std


## Calculation of uniform confidence bands ----
uniform_T_stat = lapply(1:reps, function(j) T_stat_abs[[j]]) 
uniform_T_max = lapply(1:reps, function(j) apply(uniform_T_stat[[j]], 1, max)) # max of t_stats 
uniform_q_star = lapply(1:reps, function(j) quantile(uniform_T_max[[j]], 1-alpha_sig)) # quantile of max_t_stat
uniform_CI = lapply(1:reps, function(j) list(unlist(theta_hat[[j]])-(uniform_q_star[[j]]*sigma_hat[j]),
                                             unlist(theta_hat[[j]])+(uniform_q_star[[j]]*sigma_hat[j])))
### Coverage by uniform confidence bands
coverage_true_uniform = mean(sapply(1:reps, function(k)
  sum((theta_true > uniform_CI[[k]][[1]]) & (theta_true < uniform_CI[[k]][[2]]))/nrow(X)))*100

coverage_expected_uniform = mean(sapply(1:reps, function(k)
  sum((theta_hat_expected > uniform_CI[[k]][[1]]) & (theta_hat_expected < uniform_CI[[k]][[2]]))/nrow(X)))*100

# Plotting ----
## plot with average confidence intervals -----
for (x2 in c(0.3,0.5)){
  avg_ci_lower = 0
  avg_ci_upper = 0
  avg_theta_hat = 0
  avg_T_stat_abs = 0
  avg_sigma_hat = 0
  for (i in seq(1, reps,1)){
    avg_ci_lower = avg_ci_lower + CI[[i]][[1]]
    avg_ci_upper = avg_ci_upper + CI[[i]][[2]]
    avg_theta_hat = avg_theta_hat + unlist(theta_hat[[i]])
    avg_T_stat_abs = avg_T_stat_abs + T_stat_abs[[i]]
    avg_sigma_hat = avg_sigma_hat + sigma_hat[i]
  }
  avg_ci_lower = avg_ci_lower/reps
  avg_ci_upper = avg_ci_upper/reps
  avg_theta_hat = avg_theta_hat/reps
  avg_T_stat_abs = avg_T_stat_abs/reps
  
  png(file = glue('CI_averaged_','n{formatC(as.integer(n), width=4, flag="0")}_',
                  'tau{formatC(tau*10, width=3 ,flag="0")}_',
                  'sig{formatC(sig*10, width=3 ,flag="0")}_',
                  'reps{formatC(as.integer(reps) ,flag="0")}_',
                  'x2{formatC(x2*10, width=3 ,flag="0")}_',
                  '.png'),
      width=1500, height=1500)
  par(bg='transparent')
  par(mar=c(5,6,4,1)+.1)
  pd = data.frame(X1=X$X1,X2=X$X2, theta_hat =  avg_theta_hat,
                  theta_true = theta_true, CI_L = avg_ci_lower, CI_U = avg_ci_upper)
  pd = pd[pd$X2==sprintf("%0.1f", x2),]
  plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), 
       ylim= range(pd$theta_true, pd$theta_hat, pd$CI_L, pd$CI_U),
  main= bquote(X[2]  == .(x2)),
  cex.axis = 2.5, cex.lab = 2.5, cex.main=2.5)
  lines(pd$X1, pd$theta_true ,
        ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
        col='red', main="Confidence intervals", pch=19, type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$theta_hat, col='blue', pch=19,type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
  dev.off()
}

## Plot with single/multiple confidence intervals ----
for (num in c(1,4)){
  for (x2 in c(0.3,0.5)){
    png(file = glue('CI_','n{formatC(as.integer(n), width=4, flag="0")}_',
                    'tau{formatC(tau*10, width=3 ,flag="0")}_',
                    'sig{formatC(sig*10, width=3 ,flag="0")}_',
                    'reps{formatC(as.integer(num) ,flag="0")}_',
                    'x2{formatC(x2*10, width=3 ,flag="0")}_',
                    '.png'),
        width=1500, height=1500)
    par(bg='transparent')
    par(mar=c(10,11,4,1)+1)
    par(mgp=c(8,4,0))
    plot(1, type="n", xlab=bquote(X[1]), ylab=bquote(theta), xlim=c(-0.5, 0.5),
         ylim= range(pd$theta_true, pd$theta_hat, pd$CI_L, pd$CI_U+0.1),
         main= bquote(X[2]  == .(x2)),
         cex.axis = 5, cex.lab = 5, cex.main=5)
    for (i in 1:num){
      pd = data.frame(X1=X$X1,X2=X$X2, theta_hat =  unlist(theta_hat[[i]]),
                      theta_true = theta_true, CI_L = CI[[i]][[1]], CI_U = CI[[i]][[2]] )
      pd = pd[pd$X2==sprintf("%0.1f", x2),]
      lines(pd$X1, pd$theta_true ,
            ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
            col='red', main="Confidence intervals", pch=19, type = "b", lty = 2, cex=2)
      lines(pd$X1, pd$theta_hat, col='blue', pch=19, type = "b", lty = 2, cex=2)
      lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
      lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
    }
    dev.off()
  }}

## Plot for uniform confidence bands -----
grid_T_stat = T_stat_abs[[1]] #first replication
grid_T_max = apply(grid_T_stat, 1, max) # max of t_stats 
dist_T_max = density(grid_T_max)
png(file = glue('dist',
                '.png'),
    width=1500, height=1500)

plot(dist_T_max, col= 'blue')
dev.off()

grid_q_star = quantile(grid_T_max, 1-alpha_sig) # quantile of max_t_stat
grid_CI_L = unlist(theta_hat[[1]])-(grid_q_star*sigma_hat[[1]])
grid_CI_U = unlist(theta_hat[[1]])+(grid_q_star*sigma_hat[[1]])

for (x2 in c(0.3,0.5)){
  pd = data.frame(X1=X$X1,X2=X$X2, sigma = sigma_hat[[1]],theta_hat =  unlist(theta_hat[[1]]),
                  theta_true = theta_true, CI_L = CI[[1]][[1]], CI_U = CI[[1]][[2]],
                  grid_CI_L = grid_CI_L,  grid_CI_U = grid_CI_U
  )
  #pd = pd[seq(1, nrow(pd), length.out = grids), ] #getting grids rows from our data frame
  pd = pd[pd$X2==sprintf("%0.1f", x2),]
  
  png(file = glue('CI_bands_','n{formatC(as.integer(n), width=4, flag="0")}_',
                  'tau{formatC(tau*10, width=3 ,flag="0")}_',
                  'sig{formatC(sig*10, width=3 ,flag="0")}_',
                  'grids{formatC(as.integer(grids) ,flag="0")}_',
                  'x2{formatC(x2*10, width=3 ,flag="0")}_',
                  '.png'),
      width=1500, height=1500)
  par(bg='transparent')
  par(mar=c(5,6,4,1)+.1)
  plot(1, type="n", xlab=bquote(X[1]), ylab=bquote(theta), xlim=c(-0.5, 0.5), 
       ylim=range(c(pd$theta_true,pd$theta_hat, pd$CI_L,pd$CI_U, pd$grid_CI_U,pd$grid_CI_L)),
       main= bquote(X[2]  == .(x2)),
       cex.axis = 2.5, cex.lab = 2.5, cex.main=2.5)
  lines(pd$X1, pd$theta_true ,
        ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
        col='red', main="Confidence intervals", pch=19,type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$theta_hat, col='blue', pch=19,type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$CI_L,      col='black',pch = 19,type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$CI_U,      col='black',pch = 19,type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$grid_CI_L, col='magenta',pch = 19, lty = 2, cex=1)
  lines(pd$X1, pd$grid_CI_U, col='magenta',pch = 19, lty = 2, cex=1)
  
  
  dev.off()
}
# removing undesirable variables ----
rm(current_path,ptm)

