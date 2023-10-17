library("mbsts")
library("ggplot2")
library("patchwork")
library("forecast")
library("bsts")

# Upload data
data <- readxl::read_excel("Alldata_weeklyAvg_state_2020.xlsx")

# divide 50 states into different groups
group.state <- list(
  pacific = c("Washington", "Oregon", "California", "Hawaii"),
  mountain = c("Montana", "Idaho", "Wyoming", "Nevada", "Utah", "Colorado", "Arizona", "New Mexico"),
  westnorth = c("North Dakota", "Minnesota", "South Dakota", "Nebraska", "Iowa", "Kansas", "Missouri"),
  westsouth = c("Texas", "Oklahoma", "Arkansas", "Louisiana"),
  eastsouth = c("Mississippi", "Alaska", "Tennessee", "Kentucky"),
  eastnorth = c("Illinois", "Indiana", "Wisconsin", "Michigan", "Ohio"),
  southaltantic = c("Florida", "Georgia", "South Carolina", "North Carolina", "Virginia", "West Virginia", "Maryland", "Delaware", "Alabama"),
  middlealtantic = c("Pennsylvania", "New Jersey", "New York"),
  newengland = c("Connecticut", "Rhode Island", "Massachusetts", "Vermont", "New Hampshire", "Maine")
)


# Reordering the name of state in each principal strata
for (a in 1:length(group.state)) {
  group.state[[a]] <- sort(group.state[[a]])
}


############################
# Implementing mbsts model #
############################


# Initiazing space for saving data
error.mbsts <- matrix(NA, length(group.state), 3) # Space for saving mse
tuning.parameter <- matrix(NA, 4, length(group.state)) # Optimal tuning parameter choice in mbsts model
beta.matrix.esti <- array(NA, dim = c(length(unlist(group.state)), 3, 7))
lag <- 0 # lag time
time_begin <- c(9, 21, 36) # Set up beginning time
time_end <- c(19 - lag, 34 - lag, 50 - lag) # Set up ending time
group_length <- 0

# Candidate for tuning parameter in mbsts model
rho <- 0.8
Season <- 6
vrho.space <- c(0.01, 0.05, 0.2) # candidate tuning parameter for mbsts model
lambda <- 0

# Some other parameter that will be used in training and deriving the result
alpha <- 0.05 # Type I error
threshold <- 0.5 # Feature selection, for more details, please find threshold in para.est function in mbsts package
length_X <- 51 # length of time series

# Starting training
for (i in 1:length(group.state)) {

  # data initialization
  data_subset <- subset(data, data$StateName %in% group.state[[i]]) # Subset for training
  X_init <- as.matrix(data_subset[, c(4, 6:9, 16:17)]) # covariate X in mbsts model
  number_state <- length(group.state[[i]]) # number of states in the corresponding subset
  outcome <- data_subset$caserate # Set up outcome variable
  Y <- matrix(outcome, length(outcome) / number_state) # reshape outcome
  Y <- Y * 1e-3

  # reshape covariate for mbsts model
  X <- X_init[1:length_X, ]
  array.number <- rep(NULL, number_state) # number of covariate in mbsts model
  array.number[1] <- ncol(X_init) # Initiazing array.number[1] in the above setting
  for (j in 1:(number_state - 1)) {
    X <- cbind(X, X_init[(j * length_X + 1):((j + 1) * length_X), ]) # Reshaping X for mbsts model
    array.number[j + 1] <- array.number[j] + ncol(X_init) # number of covariate
  }

  # handling missing data
  X[is.na(X) == TRUE] <- -100
  for (j in 1:nrow(X)) {
    for (k in 1:ncol(X)) {
      if (X[j, k] == -100) {
        median_value <- median(X_init[!is.na(X_init[, k %% 7]), k %% 7])
        X[j, k] <- median_value * (1 + rnorm(n = 1, sd = 1e-2)) # if X[j, k] is missing, using median
      }
    }
  }
  X <- as.matrix(X)


  # for loop for each time period
  for (k in 1:length(time_begin)) {
    Ytrain <- Y[time_begin[k]:time_end[k] + lag, ] # Training outcome, the lag effect will be considered
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

    # iteration
    number <- 0
    while (number < 10) {
      Ytrain <- Ytrain * 3 # since some range of value Y will lead to the failure of training, we time 3 for each loop to make the training valid
      number <- number + 1

      # Training
      for (s in 1:length(vrho.space)) {

        # Model initializing
        STmodel <- tsc.setting(Ytrain,
          mu = rep(1, number_state), rho = rep(rho, number_state), S = rep(Season, number_state),
          vrho = rep(vrho.space[s], number_state), lambda = rep(lambda, number_state)
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
          if (time_end[3] + lag + 1 <= 51) {
            output <- mbsts.forecast(mbsts.model,
              STmodel,
              newdata = t(as.matrix(X[time_end[k] + 1, ])),
              steps = 1
            )

            # Prediction error and its confidence interval
            error.predict <- mean(as.matrix(abs(output$pred.mean * 3^(-number) * 10^3 - Y[time_end + lag + 1, ][k, ] * 1e3))) /
              max(data$caserate[51*(0:49) + time_end[k] + lag])

            # Update optimal mbsts model and its corresponding parameter
            if (error.predict < error.optimal) {
              mbsts.model.optimal <- mbsts.model # save the updated mbsts model
              output.optimal <- output # save the updated mbsts output
              error.optimal <- error.predict # save the updated optimal error
              tuning.parameter[, i] <- c(rho, Season, vrho.space[s], lambda) # save the updated optimal hyperparameter
              beta.matrix.esti[(group_length + 1):(group_length + length(group.state[[i]])), k, ] <- matrix(rowMeans(mbsts.model.optimal@beta.hat), ncol = 7) # save the estimated beta value
            }
          }
        }
      }
    }
    error.mbsts[i, k] <- error.optimal # finding the optimal error
  }
  group_length <- group_length + length(group.state[[i]])
}


###########################
# Implementing bsts model #
###########################


# unlist state name
lag <- 0
group.state <- unlist(group.state)
data_subset <- subset(data, data$StateName %in% group.state)
error.bsts <- matrix(0, length(time_begin), length(group.state))

# bring the initial.claims data into scope
for (i in 1:length(group.state)) {
  # Get the dataset for one state
  data_raw <- subset(data_subset, data_subset$StateName %in% group.state[i])

  # Solving missing data problem in our dataset
  data_raw[is.na(data_raw) == TRUE] <- -100
  for (j in 1:nrow(data_raw)) {
    for (k in 1:ncol(data_raw)) {
      if (data_raw[j, k] == -100) {
        tmp <- unlist(data_subset[, k])
        data_raw[j, k] <- median(tmp[!is.na(tmp)]) # if data_raw[j, k] is missing, using median
      }
    }
  }

  data_raw <- data_raw[, c(4:9, 16:17)] # Choosing the covariate we will be used as regressors
  namespace <- colnames(data_raw) # Get the column name of our dataset

  # Reset the name of some regressors in case of the potential issue in formula for bsts model
  namespace[8] <- "negative"
  namespace[9] <- "positive"

  # Reinitiazing the namespace of dataset
  colnames(data_raw) <- namespace

  # bsts model training and foreasting
  for (k in 1:length(time_begin)) {

    # Initiazing the training covariate and training outcome
    data_raw_period <- data_raw[time_begin[k]:time_end[k], -2]
    data_raw_period <- cbind(data_raw[time_begin[k]:time_end[k], -2], data_raw[time_begin[k]:time_end[k] + lag, 2])

    # Initiazing bsts model
    ss <- AddLocalLinearTrend(list(), data_raw_period$caserate) # Add local linear trend component for bsts model
    ss <- AddSeasonal(ss, data_raw_period$caserate, nseasons = 6) # Add seasonal component for bsts model
    data_raw_period <- data.frame(data_raw_period) # set up our dataset as data.frame structure

    # Initiazing formula for bsts model
    namespace <- colnames(data_raw)
    formula_bsts <- as.formula(paste(
      namespace[2],
      paste(" ~ ", paste(namespace[-2], collapse = "+"))
    ))

    # Training
    bsts.model <- bsts(formula_bsts, state.specification = ss, data = data_raw_period, niter = 1000)

    # Foreasting
    pred <- predict(bsts.model, newdata = data_raw[time_end[k] + 1, -2])

    # Calculate the estimated error
    error.bsts[k, i] <- as.numeric(abs(pred$mean - data_raw[time_end[k] + lag + 1, 2])) /
      max(data$caserate[51*(0:49) + time_end[k] + lag])
  }
}



#############################
# Implementing SARIMA model #
#############################

# unlist state name
lag <- 0
group.state <- unlist(group.state)
data_subset <- subset(data, data$StateName %in% group.state)
error.SARIMA <- matrix(0, length(time_begin), length(group.state))

# Seasonal Autoregressive Integrated Moving Average (SARIMA) model
for (i in 1:length(group.state)) {
  # Get the dataset for one state
  data_raw <- subset(data_subset, data_subset$StateName %in% group.state[i])

  # Solving missing data problem in our dataset
  data_raw[is.na(data_raw) == TRUE] <- -100
  for (j in 1:nrow(data_raw)) {
    for (k in 1:ncol(data_raw)) {
      if (data_raw[j, k] == -100) {
        tmp <- unlist(data_subset[, k])
        data_raw[j, k] <- median(tmp[!is.na(tmp)]) * (1 + rnorm(n = 1, sd = 0.01)) # if data_raw[j, k] is missing, using median
      }
    }
  }

  data_raw <- data_raw[, c(4:9, 16:17)] # Choosing the covariate we will be used as regressors
  namespace <- colnames(data_raw) # Get the column name of our dataset

  for (k in 1:length(time_begin)) {
    outcome <- data_raw[time_begin[k]:(time_end[k] + 1) + lag, 2]
    covariate <- data_raw[time_begin[k]:(time_end[k] + 1), -2]
    n <- (time_end[k] + 1) - time_begin[k] + 1
    ts_data <- ts(outcome^(1/3), frequency = 6) # Set frequency to 6

    covariate_train <- as.matrix(covariate[-n, ])
    covariate_test <- as.matrix(covariate[n, ])

    # Assuming you want to use the last 100 observations as a test set
    train_data <- head(ts_data, (n - 1))
    test_data <- tail(ts_data, -(n - 1))

    # Fit SARIMA with regression for the time series
    sarima_model <- auto.arima(train_data, xreg = covariate_train)

    # Assuming you have exogenous variable values for the forecast period in 'test_data'
    forecast_values <- forecast(sarima_model, h = 6, xreg = covariate_test)
    error.SARIMA[k, i] <- abs(as.vector(forecast_values$mean)^(3) - test_data^(3)) /
      max(data$caserate[51*(0:49) + time_end[k] + lag]) 
  }
}

colMeans(error.mbsts)
rowMeans(error.bsts)
rowMeans(error.SARIMA)
