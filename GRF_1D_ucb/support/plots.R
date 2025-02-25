
# Plotting ----
## Plot of coverage with bandwidths----------------------------
png(file = glue('images/',
                'bandwidth_vs_coverage_',
                'n{formatC(as.integer(n),width=4, flag="0")}_',
                '.png'),
    width=1400, height=1400, res=300)
plot(1)
plot(node_sizes,as.numeric(unlist(coverages_widths_PW)),  type= 'l', xlab = bquote(h), ylab = 'coverage', ylim= c(70,100))
lines(node_sizes,as.numeric(unlist(coverages_widths_uniform)), col='magenta')
dev.off()

## Plot of size performance with bandwidths ------------------------
png(file = glue('images/',
                'bandwidth_vs_size_',
                'n{formatC(as.integer(n),width=4, flag="0")}_',
                '.png'),
    width=1400, height=1400, res=300)
plot(1)
plot(node_sizes,size_bandwidth_uni,  type= 'l', xlab = bquote(h), ylab = 'size', col='blue')
#lines(node_sizes,as.numeric(unlist(coverages_widths_uniform)), col='magenta')
dev.off()

## Plot of size performance with grids ----------------------------
png(file = glue('images/',
                'more_reps_grids_vs_size_pw_',
                'n{formatC(as.integer(n),width=4, flag="0")}_',
                '.png'),
    width=1400, height=1400, res=300)
plot(1)
plot(grids_x1_list,1-(as.numeric(unlist(coverages_widths_PW)/100)),  type= 'l', xlab = bquote(J), ylab = 'size'
     , col='blue')
#lines(node_sizes,as.numeric(unlist(coverages_widths_uniform)), col='magenta')
dev.off()

## Plot of Power with changing theta_0 ----------------------------
png(file = glue('images/',
                'power_curve_',
                'n{formatC(as.integer(n),width=4, flag="0")}_',
                '.png'),
    width=1400, height=1400, res=300)
plot(1)
plot(theta_0s,as.numeric(unlist(coverages_widths_PW)),  type= 'l', xlab = bquote(theta_0), ylab = 'power', ylim= c(70,100))
lines(theta_0s,as.numeric(unlist(coverages_widths_uniform)), col='magenta')
dev.off()

# Plot of density of t stat of a single point vs standard normal ----------------------------
d = density(unlist(T_stat[[1]][1,]), n=b)

fn = glue('images/',
  'test_stat_one_point',
  '.png'
)

png(file = fn)

plot(d, ylab='', xlab='',ylim = c(0,0.45), main="Density plot of bootsrap test statistic")
ylim = c(0,1)
x_std = rnorm(n, mean = 0, sd= 1)
std_d = density(x_std, n=b)
lines(std_d, col = 'blue')
dev.off()



## plot with average confidence intervals -----------------------------
# required variables (x2_fixed,CI,theta_true,theta_hat, T_stat_abs,sigma_hat,reps)
for (x2 in x2_fixed){
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
  
  png(file = glue('images/',
                  'CI_averaged_',
                  'n{formatC(as.integer(n), width=4, flag="0")}_',
                  'tau{formatC(tau*10, width=3 ,flag="0")}_',
                  'sig{formatC(sig*10, width=3 ,flag="0")}_',
                  'reps{formatC(as.integer(reps) ,flag="0")}_',
                  'x2{formatC(x2*10, width=3 ,flag="0")}_',
                  '.png'),
      width=1500, height=1500)
  
  pd = data.frame(X1=X$X1,X2=X$X2, theta_hat =  avg_theta_hat,
                  theta_true = theta_true, CI_L = avg_ci_lower, CI_U = avg_ci_upper)
  pd = pd[pd$X2==sprintf("%0.1f", x2),]
  plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), 
       ylim= range(pd$theta_true, pd$theta_hat, pd$CI_L, pd$CI_U))
  lines(pd$X1, pd$theta_true ,
        ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
        col='red', main="Confidence intervals", pch=19, type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$theta_hat, col='blue', pch=19,type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$CI_L,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$CI_U,      col='black',pch = 19,  type = "b", lty = 2, cex=0.8)
  dev.off()
}

## Plot with single/multiple confidence intervals ----------------------------
for (num in c(1,4)){
  for (x2 in c(0.3,0.5)){
    png(file = glue('images/',
                    'CI_',
                    'n{formatC(as.integer(n), width=4, flag="0")}_',
                    'tau{formatC(tau*10, width=3 ,flag="0")}_',
                    'sig{formatC(sig*10, width=3 ,flag="0")}_',
                    'reps{formatC(as.integer(num) ,flag="0")}_',
                    'x2{formatC(x2*10, width=3 ,flag="0")}_',
                    '.png'),
        width=1500, height=1500)
    
    plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5),
         ylim= range(pd$theta_true, pd$theta_hat, pd$CI_L, pd$CI_U))
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

## Plot for uniform confidence bands -----------------------------
## just for simplicity, renaming test sets as original sets
X = X_test
theta_hat = theta_hat_test
theta_true = theta_true_test
w = w_test
grid_T_stat = T_stat_abs[[1]] #first replication
grid_T_max = apply(grid_T_stat, 2, max) # max of t_stats for each bootstrap
dist_T_max = density(grid_T_max, y= 'scaled')
png(file = glue('images/',
                'dist_',
                'n{formatC(as.integer(n), width=4, flag="0")}_',
                'tau{formatC(tau*10, width=3 ,flag="0")}_',
                'sig{formatC(sig*10, width=3 ,flag="0")}_',
                'reps{formatC(as.integer(reps) ,flag="0")}_',
                '.png'),
    width=2000, height=1500)
par(bg='transparent')
plot(dist_T_max, col= 'blue', lwd = 5, cex.lab=3, cex.axis=2, 
     cex.main=3, cex.sub=3, main = 'Density plot of bootstrap test statistic')
dev.off()

grid_q_star = quantile(grid_T_max, 1-alpha_sig) # quantile of max_t_stat
grid_CI_L = unlist(theta_hat[[1]])-(grid_q_star*sigma_hat[[1]])
grid_CI_U = unlist(theta_hat[[1]])+(grid_q_star*sigma_hat[[1]])

  pd = data.frame(X1=X_test$X1, sigma = sigma_hat[,1],theta_hat =  unlist(theta_hat_test[[1]]),
                  theta_true = theta_true_test, CI_L = CI[[1]][[1]], CI_U = CI[[1]][[2]],
                  grid_CI_L = uniform_CI[[1]][[1]],  grid_CI_U = uniform_CI[[1]][[2]]
  )

  png(file = glue('images/',
                  'sin_8x',
                  'CI_bands_',
                  'n{formatC(as.integer(n), width=4, flag="0")}_',
                  'tau{formatC(tau*10, width=3 ,flag="0")}_',
                  'sig{formatC(sig*10, width=3 ,flag="0")}_',
                  'grids{formatC(as.integer(grids) ,flag="0")}_',
                  '.png'),
      width=1600, height=800)
  par(bg='transparent')
  par(mar=c(5,6,4,1)+.1)
  plot(1, type="n", xlab="X", ylab=bquote(theta), xlim=c(-0.5, 0.5), #ylim = range(c(-2.5,2.5))
       ylim=range(c(pd$theta_true,pd$theta_hat, pd$CI_L,pd$CI_U, pd$grid_CI_U,pd$grid_CI_L))
       , cex.axis = 2.5, cex.lab = 2.5
       )
  lines(pd$X1, pd$theta_true ,
        ylim=range(pd$theta_true,pd$theta_hat, pd$CI_L, pd$CI_U),
        col='red', main="Confidence intervals", pch=19,type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$theta_hat, col='blue', pch=19,type = "b", lty = 2, cex=2)
  lines(pd$X1, pd$CI_L,      col='black',pch = 19,type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$CI_U,      col='black',pch = 19,type = "b", lty = 2, cex=0.8)
  lines(pd$X1, pd$grid_CI_L, col='magenta',pch = 19, lty = 2, cex=3, lwd = 3)
  lines(pd$X1, pd$grid_CI_U, col='magenta',pch = 19, lty = 2, cex=3, lwd = 3)
  
  dev.off()
# Power Curve ----------------------------
## Calculation for power ----
h0= 0 #theta_0 under null hyp
# counting coverage of CI 
# count increased when null hypothesis not rejected because theta_0 
# under null hyp is contained in the interval

count = do.call(cbind,
                lapply(1:reps, function(j)
                  ((h0 > CI[[j]][[1]]) & (h0 < CI[[j]][[2]]))))
count_std = do.call(cbind,
                    lapply(1:reps, function(j)  
                      ((h0 > CI_std[[j]][[1]]) & (h0 < CI_std[[j]][[2]]))))
# count_uniform = do.call(cbind,
#                 lapply(1:reps, function(j)
#                   ((h0 > grid_CI_L[[j]]) & (h0 < grid_CI_U[[j]]))))


power = 1-rowSums(count)/reps
power_std = 1-rowSums(count_std)/reps
# power_uniform = 1-rowSums(count_uniform)/reps

pd = data.frame(theta_true =  theta_true_test,
                power = power , power_std = power_std)
pd = pd[order(pd$theta_true),] #ordering pd wrt theta
pd = pd[pd$theta_true != 0,] #removing theta=0
pd= pd[!duplicated(round(pd$theta_true,5)), ] #removing duplicates

## Plotting PW power curve ----------------------------
png(file = glue('images/',
                'power_curve_',
                'n{formatC(as.integer(n),width=4, flag="0")}_',
                '.png'),
    width=1400, height=1400, res=300)
plot(pd$theta_true, pd$power, type= 'l', xlab = bquote(theta), ylab = 'Power', ylim = c(0,1))
lines(pd$theta_true, pd$power_std, col='blue')
dev.off()

