library(astsa)
library(fpp3) 
library(tidyverse)
library(forecast)

# READ: Code was inspired by Hyndman fpp ed 3 section 9.9. 
# I thought his code style was nice 

# Read and prep the data
library(readr)
hydro_df <- read_csv("Data/HydroData/hyd_post.txt")
hydro_ts = ts(hydro_df[,2])

# Exploratory Data Analysis -----------------------------------------------

# look at time series plot and corresponding acf 
# Notice the seasonality
plot(hydro_df, type = "l", xlab = "time", ylab = "level")
acf(hydro_df[,2], lag = 48)

hydro_tsibble <- hydro_df %>%
  rename(time = 1, level = 2) %>%
  mutate(time = yearmonth(time)) %>%  # sets the default to 1970. For the modeling below 
  as_tsibble(index = time)

# Here is a nice version of the plot using the tsibble
hydro_tsibble %>% gg_tsdisplay() + labs(title = "Water Level Time Series")


# We take a 12 month difference and we can see the natural peaks in PACF at 12 
hydro_tsibble %>% gg_tsdisplay(difference(level, 12),
                               plot_type = "partial", lag_max = 48) + 
  labs(title = "Seasonally Differenced") + 
  ylab("")

# Data is still non-stationary so we take the first difference 
hydro_tsibble %>% gg_tsdisplay(difference(level, 12) %>% difference(1),
                               plot_type = "partial", lag_max = 60) + 
  labs(title = "Seasonal + Monthly Differenced") + 
  ylab("")

# SARIMA Model Fitting -----------------------------------------------------------

# Fit any reasonable model
fit <- hydro_tsibble %>% model(
  arima1  = ARIMA(level ~ pdq(1,1,0) + PDQ(1,1,0)),
  arima2  = ARIMA(level ~ pdq(0,1,1) + PDQ(0,1,1)),
  arima3  = ARIMA(level ~ pdq(1,1,1) + PDQ(1,1,1)),
  arima4  = ARIMA(level ~ pdq(2,1,0) + PDQ(1,1,0)),
  arima5  = ARIMA(level ~ pdq(0,1,2) + PDQ(0,1,1)),
  arima6 =  ARIMA(level ~ pdq(0,1,2) + PDQ(1,1,1)),
  arima7  = ARIMA(level ~ pdq(2,1,1) + PDQ(1,1,1)),
  arima8  = ARIMA(level ~ pdq(1,1,2) + PDQ(1,1,1)),
  arima9  = ARIMA(level ~ pdq(1,1,1) + PDQ(2,1,0)),
  arima10  = ARIMA(level ~ pdq(1,1,1) + PDQ(0,1,2)),
  arima11 = ARIMA(level ~ pdq(0,1,1) + PDQ(2,1,1)),
  arima12 = ARIMA(level ~ pdq(1,1,0) + PDQ(0,1,2)),
  arima13 = ARIMA(level ~ pdq(2,1,2) + PDQ(1,1,1)),
  arima14 = ARIMA(level ~ 1 + pdq(2, 1, 1) + PDQ(4,1,1), 
                  order_constraint = p + q + P + Q <= 10 & (constant + d + D <= 10)) # for fun 
  #auto = ARIMA(level, stepwise = FALSE, approx = FALSE)
)

fit %>%  pivot_longer(everything(), names_to = "Model name",
                      values_to = "Orders")
glance(fit) %>%  arrange(AICc) %>%  select(.model:BIC)

# Check individual residual plots
fit %>%  select(arima13) %>%  gg_tsresiduals(lag=20)
fit %>%  select(arima5) %>%  gg_tsresiduals(lag=20)
fit %>%  select(arima6) %>%  gg_tsresiduals(lag=20)
fit %>%  select(arima8) %>%  gg_tsresiduals(lag=20)
fit %>%  select(arima2) %>%  gg_tsresiduals(lag=20)

# Check normality of innovations 
augment(fit) %>% 
  filter(.model == "arima13")  %>% 
  features(.innov, ljung_box, lag = 24, dof = 6)

augment(fit) %>% 
  filter(.model == "arima5")  %>% 
  features(.innov, ljung_box, lag = 24, dof = 3)

augment(fit) %>% 
  filter(.model == "arima6")  %>% 
  features(.innov, ljung_box, lag = 24, dof = 4)

augment(fit) %>% 
  filter(.model == "arima8") %>% 
  features(.innov, ljung_box, lag = 24, dof = 5)

augment(fit) %>% 
  filter(.model == "arima2") %>% 
  features(.innov, ljung_box, lag = 24, dof = 2)

# Exponential Fitting -----------------------------------------------------

# Fit the ETS model
ets_model <- hydro_tsibble %>%
  model(ETS(level))

# View the model summary
report(ets_model)

ets_model |>
  gg_tsresiduals(lag_max = 16)

# Very bad results
augment(ets_model) |>
  features(.innov, ljung_box, lag = 16) 

# Cross Validation  -------------------------------------------------------

# I want to use the tsCV function from the forecast library so I need to make functions corresponding to my models in question 
# It takes quite a long time to run the cv values so we test only the top six models 
library(tictoc)

tic("Running the SARIMA CV")
forecast_hydro_A <- function(x, h) {forecast(Arima(x, order = c(2, 1, 2), seasonal = list(order = c(1, 1, 1), period = 12)), h=h)}
cv_errors_A <- tsCV(hydro_ts, forecast_hydro_A, h = 24)
MSE_CV_A <- sum(apply(cv_errors_A, 2, function(x) mean(x^2, na.rm = TRUE)))

forecast_hydro_B <- function(x, h) {forecast(Arima(x, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 1), period = 12)), h=h)}
cv_errors_B <- tsCV(hydro_ts, forecast_hydro_B, h = 24)
MSE_CV_B <- sum(apply(cv_errors_B, 2, function(x) mean(x^2, na.rm = TRUE)))

forecast_hydro_C <- function(x, h) {forecast(Arima(x, order = c(1, 1, 2), seasonal = list(order = c(1, 1, 1), period = 12)), h=h)}
cv_errors_C <- tsCV(hydro_ts, forecast_hydro_C, h = 24)
MSE_CV_C <- sum(apply(cv_errors_C, 2, function(x) mean(x^2, na.rm = TRUE)))

forecast_hydro_D <- function(x, h) {forecast(Arima(x, order = c(0, 1, 2), seasonal = list(order = c(1, 1, 1), period = 12)), h=h)}
cv_errors_D <- tsCV(hydro_ts, forecast_hydro_D, h = 24)
MSE_CV_D <- sum(apply(cv_errors_D, 2, function(x) mean(x^2, na.rm = TRUE)))

forecast_hydro_E <- function(x, h) {forecast(Arima(x, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12)), h=h)}
cv_errors_E <- tsCV(hydro_ts, forecast_hydro_E, h = 24)
MSE_CV_E <- sum(apply(cv_errors_E, 2, function(x) mean(x^2, na.rm = TRUE)))

forecast_hydro_F <- function(x, h) {forecast(Arima(x, order = c(1, 1, 0), seasonal = list(order = c(0, 1, 2), period = 12)), h=h)}
cv_errors_F <- tsCV(hydro_ts, forecast_hydro_F, h = 24)
MSE_CV_F <- sum(apply(cv_errors_F, 2, function(x) mean(x^2, na.rm = TRUE)))

toc()

# Final Forecast ----------------------------------------------------------
# Note: The time in this plot is made up and was set above in tsibble format

forecast(fit, h=24) %>% 
  filter(.model=='arima2') %>% 
  autoplot(hydro_tsibble) +
  labs(title = "Water Level 24 Month Forecast",
       y="Water Level")

forecast_senario_1 = forecast(fit, h=24) %>% 
  filter(.model=='arima2')

