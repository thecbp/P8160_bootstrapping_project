library(tidyverse)
library(glmnet)
library(modelr)
library(matrixcalc)
library(parallel) 
library(doParallel) 
library(foreach) 
library(iterators)

# Quick little set of functions to make creating and optimizing a Logisitc LASSO model easier

optimize = function(data, resp, lambdas, start) {
  # Parameters: 
  # data : data to train a Logistic LASSO model on
  # resp : the responses associated with the data
  # lambdas : a list of lambdas to cross-validate on
  # start : starting betas to initialize the cross-validation
  
  # returns a Logistic-LASSO model for which the test MSE is the lowest
  
  path = lambda.cv(lambdas, start, data, resp, k = 5)
  lambda.min = path[which(path$lambda == min(path$lambda)),]$lambda
  
  LL.fit = LogLASSO.CD(data, resp, start, lambda = lambda.min)
  return(list(
    finalModel = LL.fit,
    lambda.min = lambda.min,
    lambda.path = path)
    )
}

predict = function(model, newdata) {
  ### Parameters:
  # model : a LogLASSO model 
  # newdata : the data that you wish you predict on
  
  # returns an array of predicted values for the new data
  
  u = exp(as.matrix(newdata) %*% model$coefficients)
  return(u / (1 + u))
}

LogLASSO.CD = function(X, y, beta, lambda, tol = 1e-5, maxiter = 1000) {
  ### Parameters:
  # X : design matrix, does include an intercept column                                                 
  # y : response variable (should be binary)                          
  # beta : starting beta coefficients to start from                   
  # lambda : constraining parameter for LASSO penalization            
  # tol : how precise should our convergence be                       
  # maxiter : how many iterations should be performed before stopping 
  
  # return a list containing the matrix of the coefficients and the iteration matrix
  
  # Convert data into matrix for easier manipulation
  X = as.matrix(X)
  names(beta) = colnames(X) # Assign original covariate names to the betas
  
  # Iteration setup
  j = 0
  work = compile(X, y, beta)
  obj = calc.obj(beta, X, y, lambda)
  diff.obj = Inf
  path = c(iter = j, obj = obj, beta)
  beta[1] = sum(work$w * (work$z - X %*% beta)) / sum(work$w)
  while (j < maxiter && diff.obj > tol) {
    
    j = j + 1
    prev.obj = obj
    prev.beta = beta
    
    # Coordinate descent through all of the betas 
    for (k in 2:length(beta)) {
      work = compile(X, y, beta)
      z.rm.k = X[,-k] %*% beta[-k]
      val = sum(work$w * X[,k] * (work$z - z.rm.k))
      beta[k] = (1/sum(work$w * X[,k]^2)) * soft(val, lambda)  
    }
    
    # Recalculate the objective
    work = compile(X, y, beta)
    obj = calc.obj(beta, X, y, lambda)
    
    # Convergence check calculation
    diff.obj = abs(prev.obj - obj) # check difference in log-likelihood
    # diff.obj = norm(as.matrix(beta - prev.beta), "F")
    prev.obj = obj
    
    # Append it to tracking matrix
    path = rbind(path, c(iter = j, obj = obj, beta))
  } 
  
  return(list(
    path = as.tibble(path),
    coefficients = beta,
    iter = j,
    obj = obj)
  )
}

soft = function(beta, gamma) {
  ### Parameters:
  # beta : the original coefficient beta from a regression
  # gamma : the desired threshold to limit the betas at
  
  # returns a single adjusted value of the original beta
  
  return(sign(beta) * max(abs(beta) - gamma, 0))
}

calc.obj = function(betas, data, resp, lambda) {
  ### Parameters:
  # intercept : the intercept term of the betas (scalar)
  # data : the associated data for each beta in betas (n x p matrix)
  # resp : the response variable of the dataset (n x 1 array)
  # betas : all the non-intercept beta coefficients (p x 1 array)
  
  # return the log-likelihood value for a logistic model
  params = compile(data, resp, betas)
  LS = (2 * nrow(data))^(-1) * sum(params$w * (params$z - (data %*% betas))^2)
  beta.penalty = lambda * sum(abs(betas))
  return(LS + beta.penalty)
}

calc.beta.norm = function(beta1, beta2) {
  ### Parameters:
  # beta1, beta2 : beta vectors to compare
  
  # returns the Frobenius norm between two beta vectors
  return(norm(as.matrix(beta1 - beta2), "F"))
}

compile = function(data, resp, betas) {
  # Helper function to contain all the calculations for coordinate logistic regression
  p = calc.cur.p(data, betas)
  w = calc.working.weights(p)
  z = calc.working.resp(data, resp, betas)
  return(list(
    p = p,
    w = w,
    z = z
  ))
}

calc.cur.p = function(data, betas) {
  ### Parameters:
  # intercept : the intercept term of the betas (scalar)
  # data : the associated data for each beta in betas (n x p matrix)
  # betas : all the non-intercept beta coefficients (p x 1 array)
  
  # return n x 1 array of current probabilities evaluated with given betas
  
  u = data %*% betas
  return(exp(u) / (1 + exp(u)))
}

calc.working.weights = function(p) {
  ### Parameters:
  # p : the working probabilities, calculated by calc.cur.p
  
  # return n x 1 array of working weights for the data
  
  # Check for coefficient divergence, adjust for fitted probabilities 0 & 1
  close.to.1  = (1 - p) < 1e-5
  close.to.0 = p < 1e-5
  w = p * (1 - p)
  w[which(close.to.1)] = 1e-5
  w[which(close.to.0)] = 1e-5
  
  return(w)
}

calc.working.resp = function(data, resp, betas) {
  ### Parameters: 
  # intercept : the intercept term of the betas (scalar)
  # data : the associated data for each beta in betas (n x p matrix)
  # resp : the response variable of the dataset (n x 1 array)
  # betas : all the non-intercept beta coefficients (p x 1 array)
  # p : the working probabilities, calculated by calc.cur.p
  
  # return n x 1 array of working responses evaluated with given betas
  p = calc.cur.p(data, betas)
  w = calc.working.weights(p)
  return((data %*% betas) + ((resp - p) / w))
}



lambda.cv = function(lambdas, start, data, resp, k = 5) {
  ### Parameters: 
  # lambdas : the sequence of lambdas that you want to create solutions for
  # start : the starting beta coefficients
  # data : the data to estimate the coefficients from
  # resp : the response variable you're trying to predict
  
  # returns a matrix of average test MSES against a sequence of given lambdas
  
  folds = crossv_kfold(data, k = k)
  path = NULL
  i = 0
  
  for (l in lambdas) {
    # Reset the storage of the fold MSEs
    fold.mses = NULL
    i = i + 1
    for (k in 1:nrow(folds)) {
      
      # Grab the specific training indexes
      train.idx = folds[k,1][[1]][[toString(k)]]$idx
      test.idx = folds[k,2][[1]][[toString(k)]]$idx
      # Split up the data into the training and test datasets
      train.X = data[train.idx,]
      test.X = data[test.idx,]
      train.y = resp[train.idx]
      test.y = resp[test.idx]
      
      # Perform the logistic-LASSO
      fit = LogLASSO.CD(X = train.X, y = train.y, beta = start, lambda = l)
      u = exp(as.matrix(test.X) %*% fit$coefficients)
      z = u / (1 + u)
      
      # Calculate the test MSE for the fold
      fold.mse = mean((test.y - z)^2)
      fold.mses = c(fold.mses, fold.mse)
    }
    fold.se = sqrt(var(fold.mses)/5)
    
    path = rbind(path, 
                 c(lambda = l, log.lambda = log(l), 
                   avg.fold.mse = mean(fold.mses), avg.fold.se = fold.se))
    print(paste("Iteration:", i, "done"))
  }
  return(as.tibble(path))
}

lambda.cv.par = function(lambdas, start, data, resp, k = 5) {
  ### Parameters: 
  # lambdas : the sequence of lambdas that you want to create solutions for
  # start : the starting beta coefficients
  # data : the data to estimate the coefficients from
  # resp : the response variable you're trying to predict
  
  # returns a matrix of average test MSES against a sequence of given lambdas
  
  nCores = 6 # change to whatever fits your specific system
  registerDoParallel(nCores)
  
  folds = crossv_kfold(data, k = k)
  
  res = foreach(i = 1:length(lambdas), .combine = rbind) %dopar% { 
    fold.mses = NULL
    for (k in 1:nrow(folds)) {
      
      # Grab the specific training indexes
      train.idx = folds[k,1][[1]][[toString(k)]]$idx
      test.idx = folds[k,2][[1]][[toString(k)]]$idx
      
      # Split up the data into the training and test datasets
      train.X = data[train.idx,]
      test.X = data[test.idx,]
      train.y = resp[train.idx]
      test.y = resp[test.idx]
      
      # Perform the logistic-LASSO
      fit = LogLASSO.CD(X = train.X, y = train.y, beta = start, lambda = lambdas[i])
      u = exp(as.matrix(test.X) %*% fit$coefficients)
      z = u / (1 + u)
      
      # Calculate the test MSE for the fold
      fold.mse = mean((test.y - z)^2)
      fold.mses = c(fold.mses, fold.mse)
    }
    fold.se = sqrt(var(fold.mses)/5)
    c(lambda = lambdas[i], log.lambda = log(lambdas[i]), avg.fold.mse = mean(fold.mses), avg.fold.se = fold.se)
  }
  return(as.tibble(res))
}
