# mbsts_tl function is to implement mbsts-tl method that mention in paper:
#
# Description for the input in mbsts-tl function
# 
# X: A (n*K)-dimensional matrix containing all candidate predictor series for 
# each target series. K=∑ k_i is the number of all candidate predictors for 
# all target series. The first k_1 variables are the set of candidate predictors 
# for the first target series, and the next k_2 variables are the set of candidate 
# predictors for the second target series, etc.
#
# Y: A (n*m)-dimensional matrix containing multiple target series, where n is 
# the number of observations and m is the number of target series. Note that 
# usually m should not be too large.
# 
# array.number: A vector of integer values denoting the acumulated number of 
# predictors for target series. For example, if there are three target series 
# where the first has 8 predictors, the second has 6 predictors, and the third 
# has 10 predictors, then the vector is c(8,14,24).
# 
# rho.space/season.space/vrho.space/lambda.space: candidate for the hyperparameter in
# mbsts model, more details for rho, season, vrho and lambda can be found in 
# mbsts_function in R package mbsts.
#
# time_begin/time_end: the starting time and ending time for each time period.
#
# 
# Description for the output in mbsts-tl function
# 
#
# error: error is a vector for calculating normalized absolute error for each
# time period.
#
# optimal.parameter: optimal.parameter is a m*4-dimension matrix that the optimal 
# hyperparameter we use in mbsts model, where m is the number of period we explore.
# 
# beta: beta is a m*dimension matrix where dimension is the dimension of covariate 
# for all time series.
# 
# 
# Example 
# result <- mbsts_tl(X = X, array.number = array.number, Y = Y, 
#               time_begin = time_begin, time_end = time_end)
# result$error
#
#
#
library(mbsts)
mbsts_tl <- function(X, array.number, Y, lag = 0, rho.space = c(0.01, 0.2, 0.8), 
                     season.space = c(6, 10), vrho.space = c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9), 
                     lambda.space = c(0, 1), time_begin, time_end) {
  
  # Check the compatibility of covariate and outcome
  try(if (max(array.number) > nrow(X)) stop("the dimension of X does not match array.number!"))
  try(if (nrow(Y) != nrow(X)) stop("the length of X does not match the length of Y!"))
  try(if (nrow(time_begin) != nrow(time_end)) stop("the length of X does not match the length of Y!"))
  
  # Check the compatibility of grid of parameter
  try(if (min(rho.space) < 0) stop("the minimum value of rho space is less than zero!"))
  try(if (max(rho.space) > 1) stop("the maximum value of rho space is greater than one!"))
  try(if (all(season.space > 0 & season.space == round(v))) stop("the value of season space is either negative or non-integer!"))
  try(if (min(vrho.space) < 0) stop("the minimum value of vrho space is less than zero!"))
  try(if (max(vrho.space) > 1) stop("the maximum value of vrho space is greater than one!"))
  try(if (min(lambda.space) < 0) stop("the minimum value of lambda space is less than zero!"))
  try(if (max(lambda.space) > pi) stop("the maximum value of lambda space is greater than one!"))
  
  # Open the space to save the output
  time_length <- nrow(Y)
  optimal.parameter <- matrix(NA, 4, length(time_begin))
  beta.esti <- matrix(NA, ncol(X), length(time_begin))
  error.mbsts <- rep(NA, length(time_begin))
  
  # Initiazing the covariate X and Y
  X[is.na(X) == TRUE] <- -1e8
  for (j in 1:nrow(X)) {
    for (k in 1:ncol(X)) {
      if (X[j, k] == -1e8) {
        median_value <- median(X[, k])
        X[j, k] <- median_value * (1 + rnorm(n = 1, sd = 1e-2)) # if X[j, k] is missing, using median
      }
    }
  }
  X <- as.matrix(X)


  # for loop for each time period
  grid_candidate <- list(rho.space = rho.space, 
                         season.space = season.space, 
                         vrho.space = vrho.space, 
                         lambda.space = lambda.space)
  optimal_parameter <- c(grid_candidate[[1]][1], grid_candidate[[2]][1], grid_candidate[[3]][1], grid_candidate[[4]][1])

  for (k in 1:length(time_begin)) {
    Xtrain <- X[time_begin[k]:time_end[k], ] # Training covariate

    for (l in 1:dim(Xtrain)[2]) {
      if (length(unique(Xtrain[, l])) == 1) {
        Xtrain[, l] <- Xtrain[, l] * (1 + rnorm(length(Xtrain[, l]), sd = 1e-2)) # If one line is absolutely the same, add some noise in case of error in mbsts function
      }
    }
    
    # Initializing some notation to obtain optimal mbsts model
    mbsts.model.optimal <- FALSE # initializing mbsts.model.old as FALSE
    error.optimal <- 1e8 # initializing error.old as 1e8
    output.optimal <- NULL # initializing output.old as NULL
    
    for (t in 1:4) {
      Ytrain <- Y[time_begin[k]:time_end[k] + lag, ] # Training outcome, the lag effect will be considered
      scaling_parameter_Y <- max(Ytrain)
      Ytrain <- Ytrain * 1e-2 / scaling_parameter_Y

      # iteration
      number <- 0
      while (number < 10) {
        # since some range of value Y will lead to the failure of training, we time 3 for each loop to make the training valid
        Ytrain <- Ytrain * 3
        number <- number + 1

        # Choosing the optimal rho
        for (s in 1:length(grid_candidate[[t]])) {

          # Candidate parameter for mbsts model
          candidate_parameter <- rep(NA, 4)
          candidate_parameter[-t] <- optimal_parameter[-t]
          candidate_parameter[t] <- grid_candidate[[t]][s]

          # Model initializing
          STmodel <- tsc.setting(Ytrain,
            mu = rep(1, number_state),
            rho = rep(candidate_parameter[1], ncol(Ytrain)),
            S = rep(candidate_parameter[2], ncol(Ytrain)),
            vrho = rep(candidate_parameter[3], ncol(Ytrain)),
            lambda = rep(candidate_parameter[4], ncol(Ytrain))
          )

          # Training model
          mbsts.model <- tryCatch(mbsts.model <- mbsts_function(
            Y = Ytrain,
            Xtrain = Xtrain,
            STmodel = STmodel,
            ki = array.number,
            pii = matrix(rep(0.5, dim(Xtrain)[2]), nrow = dim(Xtrain)[2]),
            v0 = 50,
            mc = 700,
            burn = 200
          ), error = function(e) {
            skip_to_next <<- FALSE
          })

          # Prediction using mbsts model
          if (!isFALSE(mbsts.model)) {
            if (time_end[k] + lag + 1 <= time_length) {
              output <- mbsts.forecast(mbsts.model, STmodel,
                steps = 1,
                newdata = t(as.matrix(X[time_end[k] + 1, ]))
              ) # predict the time series using mbsts model

              # Prediction error and its confidence interval
              error.predict <- mean(as.matrix(abs(output$pred.mean * 3^(-number) / (1e-2 / scaling_parameter_Y) 
                                                  - Y[time_end[k] + lag + 1, ]))) /
                  max(Y[time_end[k] + lag, ]) # calculate the normalized absolute error

              # Update optimal mbsts model and its corresponding parameter
              if (error.predict < error.optimal) {
                mbsts.model.optimal <- mbsts.model # save the updated mbsts model
                output.optimal <- output # save the updated mbsts output
                error.optimal <- error.predict # save the updated optimal error
                optimal_parameter[t] <- grid_candidate[[t]][s] # save the updated optimal hyperparameter
              }
            }
          }
        }
      }
    }
    
    # Save the result
    beta.esti[, k] <- matrix(rowMeans(mbsts.model.optimal@beta.hat) / (1e-2 / scaling_parameter_Y), ncol = 7) # save the estimated beta value
    optimal.parameter[, k] <- optimal_parameter # Saving the optimal parameter in period k
    error.mbsts[k] <- error.optimal # finding the optimal error in period k
  }
  return(list(error = error.mbsts, optimal.parameter = optimal.parameter, beta = beta.esti))
}
