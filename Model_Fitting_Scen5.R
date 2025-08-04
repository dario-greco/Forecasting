library(astsa)
library(fpp3) 
library(tidyverse)
library(forecast)
library(fable)

# Get the pollution data
city1 = as.ts(read.csv("Data/PollutionData/pollutionCity1.txt")[,2])
city2 = as.ts(read.csv("Data/PollutionData/pollutionCity2.txt")[,2])
city3 = as.ts(read.csv("Data/PollutionData/pollutionCity3.txt")[,2])

# Had to make up some date
# Go to the tsibble world 
start_date <- ymd_hm("2024-01-01 00:00") 
time_seq <- seq(from = start_date, by = "30 mins", length.out = length(city1))

# Inspect the time series 
pollution_ts_1 <- tibble(
  datetime = time_seq,
  pollution = city1
) %>%
  as_tsibble(index = datetime)

pollution_ts_1 %>% gg_tsdisplay(pollution, plot_type = "histogram") +
  labs(title = "Pollution Level City 1")

pollution_ts_2 <- tibble(
  datetime = time_seq,
  pollution = city2
) %>%
  as_tsibble(index = datetime)

pollution_ts_2 %>% gg_tsdisplay(pollution, plot_type = "histogram") + 
  labs(title = "Pollution Level City 2")

pollution_ts_3 <- tibble(
  datetime = time_seq,
  pollution = city3
) %>%
  as_tsibble(index = datetime)

pollution_ts_3 %>% gg_tsdisplay(pollution, plot_type = "histogram") + 
  labs(title = "Pollution Level City 3")

# We can use the STL function to look at the seasonal decomposition of the time series
pollution_ts_1  %>% 
  model(
    STL(pollution ~ season(period = 48) +
          season(period = 336),
        robust = TRUE)
  ) %>% 
  components() %>% 
  autoplot() + labs(title = "STL City 1", x = "Pollution Level")

pollution_ts_2 %>% 
  model(
    STL(pollution ~ season(period = 48) +
          season(period = 336),
        robust = TRUE)
  ) %>% 
  components() %>% 
  autoplot() + labs(title = "STL City 2", x = "Pollution Level")

pollution_ts_3 %>% 
  model(
    STL(pollution ~ season(period = 48) +
          season(period = 336),
        robust = TRUE)
  ) %>% 
  components() %>% 
  autoplot() + labs(title = "STL City 3", x = "Pollution Level")

# First Time Series -------------------------------------------------------

# Get some baseline model 
fit_baseline_sarima_1 = auto.arima(pollution_ts_1)
summary(fit_baseline_sarima_1)
fit_baseline_exp_1 = pollution_ts_1 %>% model(ETS(pollution))
glance(fit_baseline_exp_1)

# We want to try some of the dynamic harmonic regression models to capture the two modes of seasonality
daily_terms <- 1:6
weekly_terms <- 1:6

param_grid <- expand.grid(K_daily = daily_terms, K_weekly = weekly_terms)
results <- param_grid %>%
  mutate(AICc = NA_real_)
for (i in 1:nrow(param_grid)) {
  
  K_daily <- param_grid$K_daily[i]
  K_weekly <- param_grid$K_weekly[i]
  
  
  fit <- pollution_ts_1 %>%
    model(DHR = ARIMA(pollution ~ fourier(48, K = K_daily) +
                        fourier(336, K = K_weekly) + trend()))
  
  model_aicc <- glance(fit) %>% pull(AICc)
  
  # Store result
  results$AICc[i] <- model_aicc
  
  cat("Fitted K_daily =", K_daily, "K_weekly =", K_weekly, "-> AICc:", model_aicc, "\n")
  
}

# Take the top three models in regards to AIC

fit_dhr1_1 <- pollution_ts_1 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 4) +
                      fourier(336, K = 2) + trend()))

fit_dhr1_2 <- pollution_ts_1 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 3) +
                      fourier(336, K = 1) + trend()))

fit_dhr1_3 <- pollution_ts_1 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 3) +
                      fourier(336, K = 2) + trend()))


# Try some SARIMA that are around the auto arima parameter s

city1_sarima_models <- pollution_ts_1 %>%
  model(
    M1 = ARIMA(pollution ~ 1 + pdq(4,0,1) + PDQ(0,0,2,48)),
    M2 = ARIMA(pollution ~ 1 + pdq(3,0,1) + PDQ(0,0,2,48)),
    M3 = ARIMA(pollution ~ 1 + pdq(4,0,2) + PDQ(0,0,1,48)),
    M4 = ARIMA(pollution ~ 1 + pdq(2,0,2) + PDQ(0,0,2,48)),
    M5 = ARIMA(pollution ~ 1 + pdq(4,0,1) + PDQ(1,0,1,48))
  )

# The best models are M1 and M2 for AICc


# CROSS VALIDATION CITY 1 -------------------------------------------------
# We do expanding window, moving the window obver 336 observations each time 
h <- 336
n_total <- nrow(pollution_ts_1)
n_windows <- 5

results <- list()

for (i in 0:(n_windows - 1)) {
  
  end_train <- n_total - h * (n_windows - i)
  start_test <- end_train + 1
  end_test <- start_test + h - 1
  
  train_ts <- pollution_ts_1[1:end_train, ]
  test_ts  <- pollution_ts_1[start_test:end_test, ]
  
  fit_models <- train_ts %>%
    model(
      auto_arima = ARIMA(pollution ~ 1 + pdq(4,0,1) + PDQ(0,0,2,48)),
      DHR1 = ARIMA(pollution ~ fourier(48, 4) + fourier(336, 2) + trend()),
      DHR2 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 1) + trend()),
      DHR3 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 1) + trend()),
      DHR4 = ARIMA(pollution ~ fourier(48, 6) + fourier(336, 4) + trend()),
      manual_arima = ARIMA(pollution ~ 1 + pdq(3,0,1) + PDQ(0,0,2,48))
    )
  
  fc <- forecast(fit_models, h = h)
  acc <- accuracy(fc, test_ts) %>% mutate(window = i + 1) # Gives back RMSE which is fine 
  
  results[[i + 1]] <- acc
  print(paste("Done window", i + 1))
}

final_results <- bind_rows(results)

summary_results <- final_results %>%
  select(-window) %>%
  group_by(.model) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Look at the results. All teh models did reasonable 
print(summary_results)

# City 2 ------------------------------------------------------------------
# We follow a similar method from above
# Get some baseline model 
fit_baseline_sarima_2 = auto.arima(pollution_ts_2)
summary(fit_baseline_sarima_2)
fit_baseline_exp_2 = pollution_ts_2 %>% model(ETS(pollution))
glance(fit_baseline_exp_2)

# We want to try some of the dynamic harmonic regression models to capture the two modes of seasonality
daily_terms <- 1:6
weekly_terms <- 1:6

param_grid <- expand.grid(K_daily = daily_terms, K_weekly = weekly_terms)
results2 <- param_grid %>%
  mutate(AICc = NA_real_)
for (i in 1:nrow(param_grid)) {
  
  K_daily <- param_grid$K_daily[i]
  K_weekly <- param_grid$K_weekly[i]
  
  
  fit <- pollution_ts_2 %>%
    model(DHR = ARIMA(pollution ~ fourier(48, K = K_daily) +
                        fourier(336, K = K_weekly) + trend()))
  
  model_aicc <- glance(fit) %>% pull(AICc)
  
  # Store result
  results2$AICc[i] <- model_aicc
  
  cat("Fitted K_daily =", K_daily, "K_weekly =", K_weekly, "-> AICc:", model_aicc, "\n")
  
}

# Get some of the best models for AIC
fit_dhr2_1 <- pollution_ts_2 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 3) +
                      fourier(336, K = 1) + trend()))

fit_dhr2_2 <- pollution_ts_2 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 4) +
                      fourier(336, K = 1) + trend()))

fit_dhr2_3 <- pollution_ts_2 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 3) +
                      fourier(336, K = 3) + trend()))


# Try some SARIMA that are around the auto arima

city2_sarima_models <- pollution_ts_2 %>%
  model(
    M1 = ARIMA(pollution ~ 1 + pdq(3,1,3) + PDQ(0,0,1,48)),
    M2 = ARIMA(pollution ~ 1 + pdq(3,1,2) + PDQ(0,0,1,48)),
    M3 = ARIMA(pollution ~ 1 + pdq(2,1,3) + PDQ(0,0,1,48)),
    M4 = ARIMA(pollution ~ 1 + pdq(2,1,2) + PDQ(0,0,1,48)),
    M5 = ARIMA(pollution ~ 1 + pdq(4,1,3) + PDQ(0,0,1,48))
  )

# Cross Validation Time Series 2 ------------------------------------------
results2 <- list()

for (i in 0:(n_windows - 1)) {
  
  end_train <- n_total - h * (n_windows - i)
  start_test <- end_train + 1
  end_test <- start_test + h - 1
  
  train_ts2 <- pollution_ts_2[1:end_train, ]
  test_ts2  <- pollution_ts_2[start_test:end_test, ]
  
  fit_models2 <- train_ts2 %>%
    model(
      auto_arima = ARIMA(pollution),
      DHR1 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 1) + trend()),
      DHR2 = ARIMA(pollution ~ fourier(48, 4) + fourier(336, 1) + trend()),
      DH32 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 3) + trend()),
      manual_arima = ARIMA(pollution ~ 1 + pdq(4,1,3) + PDQ(0,0,1,48))
    )
  
  fc2 <- forecast(fit_models2, h = h)
  acc2 <- accuracy(fc2, test_ts2) %>% mutate(window = i + 1)
  
  results2[[i + 1]] <- acc2
  print(paste("Done window", i + 1))
}

# Combine results from all windows
final_results2 <- bind_rows(results2)

# Summarize accuracy over all windows
summary_results2 <- final_results2 %>%
  select(-window) %>%
  group_by(.model) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Print average performance
print(summary_results2)

# City 3 ------------------------------------------------------------------

fit_baseline_sarima_3 = auto.arima(pollution_ts_3)
summary(fit_baseline_sarima_3)
fit_baseline_exp_3 = pollution_ts_3 %>% model(ETS(pollution))
glance(fit_baseline_exp_3)

# We want to try some of the dynamic harmonic regression models to capture the two modes of seasonality
daily_terms <- 1:6
weekly_terms <- 1:6

param_grid <- expand.grid(K_daily = daily_terms, K_weekly = weekly_terms)
results3 <- param_grid %>%
  mutate(AICc = NA_real_)
for (i in 1:nrow(param_grid)) {
  
  K_daily <- param_grid$K_daily[i]
  K_weekly <- param_grid$K_weekly[i]
  
  
  fit <- pollution_ts_3 %>%
    model(DHR = ARIMA(pollution ~ fourier(48, K = K_daily) +
                        fourier(336, K = K_weekly)))
  
  model_aicc <- glance(fit) %>% pull(AICc)
  
  # Store result
  results3$AICc[i] <- model_aicc
  
  cat("Fitted K_daily =", K_daily, "K_weekly =", K_weekly, "-> AICc:", model_aicc, "\n")
  
}
# These were the top modles with respect to AIC
fit_dhr3_1 <- pollution_ts_3 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 3) +
                      fourier(336, K = 1) ))

fit_dhr3_2 <- pollution_ts_3 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 5) +
                      fourier(336, K = 1)))


fit_dhr3_3 <- pollution_ts_3 %>%
  model(DHR = ARIMA(pollution ~ fourier(48, K = 3) +
                      fourier(336, K = 2)))

# Try some SARIMA that are around teh auto arima parameter s

city3_sarima_models <- pollution_ts_3 %>%
  model(
    M1 = ARIMA(pollution ~ 1 + pdq(1,1,0) + PDQ(0,0,2,48)),
    M2 = ARIMA(pollution ~ 1 + pdq(1,1,1) + PDQ(0,0,2,48)),
    M3 = ARIMA(pollution ~ 1 + pdq(1,1,1) + PDQ(0,0,2,48)),
    M4 = ARIMA(pollution ~ 1 + pdq(2,1,0) + PDQ(0,0,2,48)),
    M5 = ARIMA(pollution ~ 1 + pdq(2,1,1) + PDQ(0,0,2,48))
  )
city3_sarima_models %>% glance()
# Cross Validation Time Series 3 ------------------------------------------
h <- 336
n_total <- nrow(pollution_ts_3)
n_windows <- 5

results3 <- list()

for (i in 0:(n_windows - 1)) {
  
  end_train <- n_total - h * (n_windows - i)
  start_test <- end_train + 1
  end_test <- start_test + h - 1
  
  train_ts3 <- pollution_ts_3[1:end_train, ]
  test_ts3  <- pollution_ts_3[start_test:end_test, ]
  
  fit_models3 <- train_ts3 %>%
    model(
      auto_arima = ARIMA(pollution ~ 1+  pdq(1,1,0) + PDQ(0,0,2,48)),
      DHR1 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 2)),
      DHR2 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 1)),
      DHR3 = ARIMA(pollution ~ fourier(48, 5) + fourier(336, 2)),
      manual_arima = ARIMA(pollution ~ 1 + pdq(1,1,1) + PDQ(0,0,2,48))
    )
  
  fc3 <- forecast(fit_models3, h = h)
  acc3 <- accuracy(fc3, test_ts3) %>% mutate(window = i + 1)
  
  results3[[i + 1]] <- acc3
  print(paste("Done window", i + 1))
}


final_results3 <- bind_rows(results3)

summary_results3 <- final_results3 %>%
  select(-window) %>%
  group_by(.model) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")


print(summary_results3)

# FORECASTS ---------------------------------------------------------------
# City 1
best_model_1 <- pollution_ts_1 %>%
  model(DHR2 = ARIMA(pollution ~ fourier(48, 3) + fourier(336, 1) + trend()))

fc_best <- forecast(best_model_1, h = 336)

fc_best %>%
  autoplot(pollution_ts_1) +
  ggtitle("Forecast from Best Model (DHR2)") +
  xlab("Time") + ylab("Pollution Level")

# City 2
best_model_2 <- pollution_ts_2 %>%
  model(auto_arima = ARIMA(pollution))

fc_best2 <- forecast(best_model_2, h = 336)

fc_best2 %>%
  autoplot(pollution_ts_2) +
  ggtitle("Forecast from Best Model ARIMA AUTO") +
  xlab("Time") + ylab("Pollution Level")

# City 3
best_model_3 <- pollution_ts_3 %>%
  model(DHR2 = ARIMA(pollution ~ fourier(48, 3) +
                       fourier(336, 1)))

fc_best3 <- forecast(best_model_3, h = 336)

fc_best3 %>%
  autoplot(pollution_ts_3) +
  ggtitle("Forecast from Best Model (DHR2)") +
  xlab("Time") + ylab("Pollution Level")


