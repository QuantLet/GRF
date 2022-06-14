[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **GRF_mindscholars_causalRF** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: GRF_mindscholars_causalRF

Published in: 'METIS'

Description: 'Computation and visualisation of (mean) ICE (individual conditional expectation) with and without treatment effect computed with a generalized Random Forest, for variable S3 (=Students self-reported expectation for success in the future) in the National Mindset Study data set used by Athey and Wager in "Estimating Treatment Effects with Causal Forests: An Application".'

Keywords: 'causal, treatment effect, GRF, generalized random forest, ICE, individual conditional expectation, approximation, mindset, mindset nudging'

Author: 'Marius Sterling'

See also: ''

Submitted:  '31.08.2020'

```

![Picture1](ICE_mean.png)

![Picture2](ICE_one.png)

### R Code
```r

# Example based on 
# https://github.com/grf-labs/grf/tree/master/experiments/

set.seed(42)
rm(list = ls())
libraries = c('grf', 'pROC', 'caret', 'future.apply', 'Hmisc',
              'lmtest', 'ggplot2', 'glue')
lapply(libraries,function(x)if(!(x %in% installed.packages())){
  install.packages(x)})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Defining sink -----------------------------------------------------------
sink_on = function() sink(file = 'log.txt', append = TRUE)
sink_off = function() sink(file = NULL)

# sink_on = function() return()
# sink_off = function() return()


# Defining rsq and adj rsq ------------------------------------------------

rsq = function(x, x_pred){
  ss_tot = sum((x - mean(x))**2)
  ss_res = sum((x - x_pred)**2)
  return(1- ss_res/ss_tot)
}
rsq_adj = function(x, x_pred, p, model_with_const=TRUE){
  n = length(x)
  return(1 - (1-rsq(x, x_pred)) * (n - 1) / (n - p - 1 - ifelse(model_with_const, 1, 0)))
}

# Read in data ------------------------------------------------------------
data.all = read.csv("synthetic_data.csv")
data.all$schoolid = factor(data.all$schoolid)

DF = data.all[,-1]
school.id = as.numeric(data.all$schoolid)
school.mat = model.matrix(~ schoolid + 0, data = data.all)
school.size = colSums(school.mat)

# It appears that school ID does not affect pscore. So it can be ignored
# in modeling.
w.lm = glm(Z ~ ., data = data.all[,-3], family = binomial)
summary(w.lm)

# Setting up data set -----------------------------------------------------
trainIndex = createDataPartition(DF$Z, p = 0.75, list = FALSE)

W = DF$Z[trainIndex]
W_test = DF$Z[-trainIndex]
Y = DF$Y[trainIndex]
Y_test = DF$Y[-trainIndex]
X.raw = DF[,-(1:2)]

# Creation of dummy variables for multi label categorical variable C1 (student 
# race/ethnicity, 15 labels) and XC (School-lebel categorical variable for 
# urbanicity of the school)
C1_ = factor(X.raw$C1)
C1.exp = model.matrix(~ C1_ + 0)
XC_ = factor(X.raw$XC)
XC.exp = model.matrix(~ XC_ + 0)
X = cbind(X.raw[,-which(names(X.raw) %in% c("C1", "XC"))], C1.exp, XC.exp)
X_test = X[-trainIndex,]
X = X[trainIndex,]

# OLS regression ----------------------------------------------------------
model <- lm(Y ~ ., data = cbind(X, W))
Y_hat = predict(model)

sink_on()
print('linear regression')
summary(model)
print(glue::glue("MSE: {mean((Y_hat - Y)**2)}"))
print(glue::glue("rsq: {round(rsq(Y,Y_hat),4)}"))
print(glue::glue("rsq_adj: {round(rsq_adj(Y,Y_hat, model$rank),4)}"))
sink_off()

Y_test_hat = predict(model, newdata = cbind(X_test, W=W_test))
sink_on()
print('Test:')
print(glue::glue("MSE: {mean((Y_test_hat - Y_test)**2)}"))
print(glue::glue("rsq: {round(rsq(Y_test, Y_test_hat),4)}"))
print(glue::glue("rsq: {round(rsq_adj(Y_test, Y_test_hat, model$rank),4)}"))
sink_off()

# causal forest -------------------------------------------------------
sink_on()
print('')
print('regression RF:')
regRF <- regression_forest(X = X, Y = Y, tune.parameters = 'all',
                           num.trees = 2000)
save(regRF, file = 'regression_forest.Rdata')
rf = regRF
print(rf)
x = variable_importance(rf, max.depth = 50)
row.names(x) = colnames(X)
print(round(x[order(x, decreasing = TRUE),],4))

# print('propensity score RF:')
# propensityRF <- regression_forest(X = X, Y = W, tune.parameters = 'all', num.trees = 2000)
# roc = pROC::roc(W, propensityRF$predictions[,1])
# plot(roc)

print('causal RF:')
cRF <- causal_forest(X = X, W = W, Y = Y, Y.hat = regRF$predictions,
                     tune.parameters = 'all', num.trees = 2000)
save(cRF, file = 'causal_forest.Rdata')
rf = cRF
print(rf)
x = variable_importance(rf, max.depth = 50)
row.names(x) = colnames(X)
print(round(x[order(x, decreasing = TRUE),],4))



# MSE, rsq, adj-rsq on train and test -------------------------------------
Y_hat = regRF$predictions
print("train:")
print(glue::glue("MSE: {mean((Y_hat - Y)**2)}"))
print(glue::glue("rsq: {round(rsq(Y, Y_hat),4)}"))
print(glue::glue("rsq_adj: {round(rsq_adj(Y, Y_hat, ncol(regRF$X.orig)),4)}"))

Y_test_hat = predict(regRF, newdata = X_test)
print("test:")
print(glue::glue("MSE: {mean(((Y_test_hat - Y_test)**2)[,1])}"))
print(glue::glue("rsq: {round(rsq(Y_test, Y_test_hat),4)}"))
print(glue::glue("rsq_adj: {round(rsq_adj(Y_test, Y_test_hat, ncol(regRF$X.orig)),4)}"))
sink_off()


# get_Y_hat function ------------------------------------------------------
get_var_xs = function(X, var, l = 100){
  return(seq(min(X[,var]), max(X[,var]), length.out = l))
}
get_Y_hat_for_var = function(regRF = regRF, var, i = NULL,
                             comp_variance = TRUE, X = X, l = 100,
                             cRF = NULL){
  x = get_var_xs(X = X, var = var, l=l)
  X2 = data.frame(x)
  # X2 = data.frame(seq(-1, 1, length.out = 100))
  colnames(X2) = var
  if(is.null(i)){
    # X_fill = apply((X[,!colnames(X) %in% var]), 2, mean)
    X_fill = apply((X[,!colnames(X) %in% var]), 2, median)
  }else{
    X_fill = X[i, !colnames(X) %in% var]  
  }
  X2[, 2:(length(X_fill) + 1)] = X_fill
  colnames(X2)[2:ncol(X2)] = names(X_fill)
  
  Y_hat = predict(regRF, X2, estimate.variance = comp_variance)
  if(!is.null(cRF)){
    Y_hat[,'treatment_effect'] = predict(cRF, X2)
  }
  return(Y_hat)
}

# ICE (Individual Conditional Expectation) --------------------------------
var = 'S3'
Y_hat = get_Y_hat_for_var(
  regRF = regRF, var = var, i = 2,
  comp_variance = TRUE, X = X_test, l = 7,
  cRF = cRF
)
save(Y_hat, file = 'ICE_one.Rdata')
png(file='ICE_one.png', bg = 'transparent')
  x = get_var_xs(X = X_test, var = var, l = 7)
  sigma.hat = sqrt(Y_hat$variance.estimates)
  ci_l = Y_hat$predictions - 1.96 * sigma.hat
  ci_u = Y_hat$predictions + 1.96 * sigma.hat
  plot(x, Y_hat$predictions,
       xlab = var, ylab = "ICE", type = "l", col = 'blue', lwd = 2,
       # ylim = range(Y_hat$predictions, 0, 1),
       ylim = c(min(ci_l), max(ci_u + max(0, Y_hat$treatment_effect))))
  lines(x, ci_u, lty = 2, col ='blue')
  lines(x, ci_l, lty = 2, col ='blue')
  lines(x, Y_hat$predictions + Y_hat$treatment_effect, lwd = 2, col = 'red')
  lines(x, ci_u + Y_hat$treatment_effect, lty = 2, col ='red')
  lines(x, ci_l + Y_hat$treatment_effect, lty = 2, col ='red')
dev.off()
# Mean ICE ----------------------------------------------------------------

plan(multisession) ## Run in parallel on local computer
nrow(X_test)
Y_hats = future_lapply(1:nrow(X_test), function(i) {
  get_Y_hat_for_var(
    regRF = regRF,
    var = 'S3',
    i = i,
    X = X_test,
    l = 7,
    comp_variance = TRUE,
    cRF = cRF
  )},
  future.packages	= c('grf')
)
save(Y_hats, file = 'ICE_mean_test.Rdata')

Y_hat_predictions = list()
Y_hat = list()
for(w in unique(W_test)){
  Y_hat[[as.character(w)]] = data.frame(
    predictions = apply(sapply(Y_hats[W_test == w], function(x) x$predictions + w * x$treatment_effect), 1, mean),
    variance.estimates = apply(sapply(Y_hats[W_test == w], function(x) x$variance.estimates), 1, mean)
  )
  Y_hat_predictions[[as.character(w)]] = sapply(Y_hats[W_test == w], function(x) x$predictions + w * x$treatment_effect)
}

png(file='ICE_mean.png', bg = "transparent")
  var = 'S3'
  x = get_var_xs(X = X_test, var = var, l = 7)
  
  u0 = Y_hat[['0']]$predictions + 1.96 * sqrt(Y_hat[['0']]$variance.estimates)
  l0 = Y_hat[['0']]$predictions - 1.96 * sqrt(Y_hat[['0']]$variance.estimates)
  u1 = Y_hat[['1']]$predictions + 1.96 * sqrt(Y_hat[['1']]$variance.estimates)
  l1 = Y_hat[['1']]$predictions - 1.96 * sqrt(Y_hat[['1']]$variance.estimates)
  
  plot(x, Y_hat[['0']]$predictions, ylim = range(u0, u1, l0, l1, 0, 1),
       xlab = var, ylab = "Mean ICE", type = "l", col = 'blue', lwd = 2)
  lines(x, u0, lty = 2, col ='blue')
  lines(x, l0, lty = 2, col ='blue')
  
  lines(x, Y_hat[['1']]$predictions, col ='red')
  lines(x, u1, lty = 2, col ='red')
  lines(x, l1, lty = 2, col ='red')
dev.off()

```

automatically created on 2020-09-14