library("mbsts")
library("ggplot2")


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


# Sorting the name of state for each region
for (a in 1:length(group.state)) {
  group.state[[a]] <- sort(group.state[[a]])
}


# Set up for hyperparameter in MBSTS model
lag <- 0 # lag effect
length_X <- 51 # Set up the length of time
time_begin <- c(9, 21, 36) # Set up beginning time
time_end <- c(19 - lag, 34 - lag, 50 - lag) # Set up ending time

# Set up space for saving result
alpha <- 0.05
threshold <- 0.5

# data initialization
i <- 1
data_subset <- subset(data, data$StateName %in% group.state[[i]]) # Subset for training
X_init <- as.matrix(data_subset[, c(4, 6:9, 16:17)]) # covariate X in mbsts model
number_state <- length(group.state[[i]]) # number of states in the corresponding subset
outcome <- data_subset$caserate # Set up outcome variable
Y <- matrix(outcome, length(outcome) / number_state) # reshape outcome

# reshape covariate for mbsts model
X <- scale(X_init[1:length_X, ]) # Standardized covariate X
array.number <- rep(NULL, number_state) # number of covariate in mbsts model
array.number[1] <- ncol(X_init) # Initiazing array.number[1] in the above setting
for (j in 1:(number_state - 1)) {
  X <- cbind(X, X_init[(j * length_X + 1):((j + 1) * length_X), ]) # Reshaping X for mbsts model
  array.number[j + 1] <- array.number[j] + ncol(X_init) # number of covariate
}

# handling missing data
X[is.na(X) == TRUE] <- -1e4
for (j in 1:nrow(X)) {
  for (k in 1:ncol(X)) {
    if (X[j, k] == -1e4) {
      X[j, k] <- median(X[, k])(1 + rnorm(n = 1, sd = 1e-4)) # if X[j, k] is missing, using median
    }
  }
}
X <- as.matrix(X)
X <- scale(X)

# tuning parameter in mbsts model
rho <- 0.8
season <- 6
vrho <- 0.01
season <- 6
lambda <- 0

# Start training
# for loop for each time period
set.seed(12644)
for (k in 1:length(time_begin)) {
  Ytrain <- Y[time_begin[k]:time_end[k] + lag, ] # Training outcome, the lag effect will be considered
  Xtrain <- X[time_begin[k]:time_end[k], ] # Training covariate

  # Model initializing
  STmodel <- tsc.setting(Ytrain,
    mu = rep(1, number_state), rho = rep(rho, number_state), S = rep(season, number_state),
    vrho = rep(vrho, number_state), lambda = rep(lambda, number_state)
  )

  # Training model
  mbsts.model <- mbsts_function(
    Y = Ytrain,
    Xtrain = Xtrain,
    STmodel = STmodel,
    ki = array.number,
    pii = matrix(rep(0.5, dim(Xtrain)[2]), nrow = dim(Xtrain)[2]),
    v0 = 50,
    mc = 700,
    burn = 200
  )

  # Prediction using mbsts model
  assign(paste("mbsts.model", k, sep = ""), mbsts.model)
  if (time_end[3] + lag + 1 <= 51) {
    output <- mbsts.forecast(mbsts.model,
      STmodel,
      newdata = t(as.matrix(X[time_end[k] + 1, ])),
      steps = 1
    )
    pred_matrix[k, ] <- output$pred.mean
    error.matrix[k, ] <- abs(pred_matrix[k, ] - Y[time_end + lag + 1, ][k, ])
  }
}

# Plot result for mbsts model
pred_plot_data <- data.frame(
  time = 1:3,
  pred = pred_matrix[, k],
  true = Y[time_end + lag + 1, ][, k]
)
for (k in 1:number_state) {
  p <- ggplot2::ggplot(data = pred_plot_data, aes(x = time, y = pred)) +
    geom_line(aes(color = "prediction")) +
    geom_point(aes(color = "prediction")) +
    geom_line(aes(y = true, color = "true value")) +
    geom_point(aes(y = true, color = "true value")) +
    ylab(group.state[[i]][k]) +
    theme_classic()
  assign(paste("p", k, sep = ""), p)
}
grid.arrange(p1, p2, p3, p4, ncol = 2)

# parameter estimation and feature selection in mbsts model
para.est(object = mbsts.model1, prob.threshold = 0.4)$index
para.est(object = mbsts.model2, prob.threshold = 0.4)
para.est(object = mbsts.model3, prob.threshold = 0.4)
plot_prob(object = mbsts.model1, prob.threshold = 0.4, varnames = colnames(Xtrain))
