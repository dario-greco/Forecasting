library(rugarch)
library(tidyverse)
library(sn) # for skew normal error dist 
library(astsa)

# Load in Data 
fold_path <- "Data/Stock/"
stock_list_ts <- vector("list", 40)

# Loop over stock1.txt to stock40.txt
for (i in 1:40) {
  file_name <- paste0(fold_path, "stock", i, ".txt")
   stock <- read.csv(file_name)
   ts_data <- ts(stock[,2])
  stock_list_ts[[i]] <- ts_data
}

# Plot all the stocks 
for(i in 1:40){
  tsplot(stock_list_ts[[i]])
}

# All Stocks --------------------------------------------------------------
# stock_list_ts is a list of all the stocks as a time series object
best_models_vec <- vector("character", length(stock_list_ts))

# Trying to minimize this error
calc_ERR <- function(actual, forecast, alpha = 0.15) {
  mean((actual - forecast) * (alpha - (actual <= forecast)))
}

for(k in 1:40){
 # Pull one stock to test method
  one_stock <- stock_list_ts[[k]]
  # The models we will use. Standard Garch (1,1) with different error distributions. 
  garch_spec_gauss <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model = list(armaOrder = c(0,0)),
    distribution.model = "norm"
  )
  
  garch_spec_snorm <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model = list(armaOrder = c(0,0)),
    distribution.model = "snorm"
  )

  # We will fit an expanding window approach starting with 120 observations
  window_exp = seq(120, 140, 1)
  n.ahead = 10
  alpha = 0.15
  
  niter = length(window_exp)
  
  ERR_gauss <- numeric()
  ERR_gauss_np <- numeric()
  ERR_snorm_np <- numeric()
  for(i in 1:niter){
    # define our expanding window 
    train_size <- window_exp[i]
    train_data <- one_stock[1:train_size] # expanding window approach
    test_data <- one_stock[(train_size + 1):(train_size + 10)] # We want to test on ten steps ahead
    
    # Fit some models 
    garch_fit_gauss <- ugarchfit(spec = garch_spec_gauss, data = train_data)
    garch_fit_snorm <- ugarchfit(spec = garch_spec_snorm, data = train_data)
    
    # Forecast
    garch_fc_gauss <- ugarchforecast(garch_fit_gauss, n.ahead = 10)
    garch_fc_snorm <- ugarchforecast(garch_fit_snorm, n.ahead = 10)
    
    # Get the estimated vol 
    vol_fc_gauss <- sigma(garch_fc_gauss)
    vol_fc_snorm <- sigma(garch_fc_snorm)
    
    # VaR approaches
    # standard parametric for gauss
    VaR_Gauss_para <- qnorm(alpha)*vol_fc_gauss
    
    # non-parametric gauss
    gauss_resid <- residuals(garch_fit_gauss, standardize = TRUE)
    VaR_Gauss_non_para <- quantile(gauss_resid, probs = 0.15)*vol_fc_gauss
    
    # non-parametric snorm
    snorm_resid <- residuals(garch_fit_snorm, standardize = TRUE)
    VaR_snorm_non_para <- quantile(snorm_resid, probs = 0.15)*vol_fc_snorm
    
    # Get the ERR
    ERR_gauss[i] <- calc_ERR(test_data, VaR_Gauss_para)
    ERR_gauss_np[i] <- calc_ERR(test_data, VaR_Gauss_non_para)
    ERR_snorm_np[i] <- calc_ERR(test_data,VaR_snorm_non_para )
    
  }
  
  mean_ERR_gauss     <- mean(ERR_gauss, na.rm = TRUE)
  mean_ERR_gauss_np  <- mean(ERR_gauss_np, na.rm = TRUE)
  mean_ERR_snorm_np  <- mean(ERR_snorm_np, na.rm = TRUE)
  
  # Combine into a named vector
  ERR_summary <- c(
    "GARCH(1,1) Gaussian Parametric" = mean_ERR_gauss,
    "GARCH(1,1) Gaussian Nonparametric" = mean_ERR_gauss_np,
    "GARCH(1,1) Skew-Norm Nonparametric" = mean_ERR_snorm_np
  )
  
  # Find the best model (lowest ERR)
  best_model <- names(ERR_summary)[which.min(ERR_summary)]
  cat("Best model based on lowest ERR:", best_model)
  best_models_vec[k] <- best_model
  
}

# GET FORECASTS -----------------------------------------------------------
# we want to take our best model and train it on the full dataset for the 10 ahead
VaR_forecasts <- list()

# Loop through each the "best model" stock list and fit the best model accoridngly. 
for(k in 1:40) {
  
  one_stock <- stock_list_ts[[k]]
  best_model <- best_models_vec[k]
  if (best_model == "GARCH(1,1) Gaussian Parametric") {
    
    garch_fit <- ugarchfit(spec = garch_spec_gauss, data = one_stock)
    garch_fc <- ugarchforecast(garch_fit, n.ahead = 10)
    
    sigma_fc <- sigma(garch_fc)
    VaR_forecast <- qnorm(alpha) * sigma_fc
    
  } else if (best_model == "GARCH(1,1) Gaussian Nonparametric") {
    
    garch_fit <- ugarchfit(spec = garch_spec_gauss, data = one_stock)
    garch_fc <- ugarchforecast(garch_fit, n.ahead = 10)
    
    sigma_fc <- sigma(garch_fc)
    
    gauss_resid <- residuals(garch_fit, standardize = TRUE)
    q_resid <- quantile(gauss_resid, probs = alpha)
    
    VaR_forecast <- q_resid * sigma_fc
    
  } else if (best_model == "GARCH(1,1) Skew-Norm Nonparametric") {
    
    garch_fit <- ugarchfit(spec = garch_spec_snorm, data = one_stock)
    garch_fc <- ugarchforecast(garch_fit, n.ahead = 10)
    
    sigma_fc <- sigma(garch_fc)
    
    snorm_resid <- residuals(garch_fit, standardize = TRUE)
    q_resid <- quantile(snorm_resid, probs = alpha)
    
    VaR_forecast <- q_resid * sigma_fc
  }
  
  # Store forecast for stock k
  VaR_forecasts[[k]] <- VaR_forecast
  
}

# Forecast Plotting + Other Viz ---------------------------------------------------------------

# Stock 1 -----------------------------------------------------------------
# get stock 1 fc
VaR_for_stock1 <- as.vector(VaR_forecasts[[1]])
time1 <- seq(151, 160, 1)

# Add the historical vol 
garch_fit_1 <- ugarchfit(spec = garch_spec_gauss, data = stock_list_ts[[1]][1:150])
sigma_estimates_1 <- sigma(garch_fit_1)
historical_VaR_1 <- qnorm(alpha) * sigma_estimates_1  # parametric quantile
VaR_series_1 <- c(as.numeric(historical_VaR_1), as.numeric(VaR_forecasts[[1]]))


# Plot Stock 1 Vol forecast
plot(stock_list_ts[[1]], type = "l", col = "black", lwd = 1.5,
     xlab = "Time", ylab = "Log Daily Returns", main = "Stock 1: VaR Forecast",
     xlim = c(0, 160))

lines(1:160, VaR_series_1, col = "red", lwd = 2)
legend("topright", legend = c("VaR"),
       col = c( "red"), lwd = 2, bty = "n")


# Stock 2 -----------------------------------------------------------------

VaR_for_stock2 <- as.vector(VaR_forecasts[[2]])
time2 <- seq(151, 160, 1)

# plot it 
plot(stock_list_ts[[2]], ylab = "Log Daily Returns Stock 2", 
     xlab = "Time", 
     main = "VaR(0.15) Forecast for Stock 2", xlim = c(0, 160))
lines(time2, VaR_for_stock2, col = "red")

# Add the historical vol 
garch_fit_2 <- ugarchfit(spec = garch_spec_gauss, data = stock_list_ts[[2]][1:150])
sigma_estimates_2 <- sigma(garch_fit_2)
resid_2 <- residuals(garch_fit_2, standardize = TRUE)
q_resid_2 <- quantile(resid_2, probs = alpha)
historical_VaR_2 <- q_resid_2 * sigma_estimates_2  


VaR_series_2 <- c(as.numeric(historical_VaR_2), as.numeric(VaR_forecasts[[2]]))

plot(stock_list_ts[[2]], type = "l", col = "black", lwd = 1.5,
     xlab = "Time", ylab = "Log Daily Returns", main = "Stock 2: VaR Forecast",
     xlim = c(0, 160))

lines(1:160, VaR_series_2, col = "red", lwd = 2)
legend("topright", legend = c("VaR "),
       col = c("red"), lwd = 2, bty = "n")

