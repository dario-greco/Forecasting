library(astsa)
library(fpp3)  
library(tidyverse)
library(forecast)
library(imputeTS)

# Scenario 4 ---------------------------------------------------------------

# Read in the data 
prod_1 = read.csv("Data/ProdData/prod_1.txt")
prod_1 = prod_1 %>% 
  select(V1, V2) %>% 
  rename(time = V1, Car = V2) %>% 
  mutate(time = yearmonth(time)) %>% 
  as_tsibble(index = time)

prod_1 %>% gg_tsdisplay(plot_type = c("auto"))
prod_1 %>% autoplot(Car) + labs(title = "Car Production")

prod_2 = read.csv("Data/ProdData/prod_2.txt")
prod_2 = prod_2 %>% 
  select(V1, V2) %>% 
  rename(time = V1, Steel = V2) %>% 
  mutate(time = yearmonth(time)) %>% 
  as_tsibble(index = time)

prod_2 %>% autoplot(Steel) + labs(title = 'Steel Production')
prod_2 %>% gg_tsdisplay(log(Steel),plot_type = c("auto"))

eng_1 = read.csv("Data/ProdData/eng_1.txt")
eng_1 = eng_1 %>% 
  select(V1, V2) %>% 
  rename(time = V1, Gas = V2) %>% 
  mutate(time = yearmonth(time)) %>% 
  as_tsibble(index = time)

eng_1 %>% gg_tsdisplay(Gas, (plot_type = c("auto")))
eng_1 %>% autoplot(Gas) + labs(title = "Gas Consumption")
eng_1 %>% gg_tsdisplay(log(Gas), (plot_type = c("auto")))

eng_2 = read.csv("Data/ProdData/eng_2.txt")
eng_2 = eng_2 %>% 
  select(V1, V2) %>% 
  rename(time = V1, Elec = V2) %>% 
  mutate(time = yearmonth(time)) %>% 
  as_tsibble(index = time)


eng_2 %>% gg_tsdisplay(Elec,  plot_type = "auto")
eng_2 %>% autoplot(Elec) + labs(title = "Electricity Consumption")
eng_2 %>% gg_tsdisplay(log(Elec),  plot_type = "auto")

# And let us recall our orginal time series
ggplot_na_imputations(beer_ts, fit2_imp)
autoplot(log(fit2_imp))

# Model Fitting -----------------------------------------------------------
# With all the time series data in these scenario we see an increase in variance 
# We want to adjust for that by taking the log. The rest of the non-stationarity will be 
# dealt with by taking the difference. 

# First we fit some regression models to see what the best predictors could be 

# Regression Fitting ------------------------------------------------------

y = log(fit2_imp)

Elec = log(eng_2$Elec)    
Gas  = log(eng_1$Gas)    
Steel = log(prod_2$Steel) 

xreg1 = cbind(
  Elec  = Elec,
  Gas = Gas,
  Steel = Steel
)

xreg2 = cbind(Gas = Gas,
              Steel = Steel)
xreg3 = cbind(Elec = Elec,
              Elec2 = Elec^2)
xreg4 = cbind(Steel = Steel,
              Elec = Elec)
xreg5 = cbind(Elec = Elec, 
              Gas = Gas)

model_lm1 = lm(y ~ xreg[,"Elec"] + I(xreg[, "Elec"]^2))
summary(model_lm1)

model_lm2 = lm(y ~ xreg[, "Gas"])
summary(model_lm2)

model_lm3 = lm(y ~ xreg[, "Steel"])
summary(model_lm3)

model_lm4 = lm(y ~ xreg)
summary(model_lm4)

model_lm5 = lm(y ~ xreg[,"Gas"] + xreg[, "Steel"])
summary(model_lm5)

log_Steel = as.vector(log(prod_2$Steel))
log_Gas = as.vector(log(eng_1$Gas))
log_Elec = as.vector(log(eng_2$Elec))

# Insepct the residuals from some of the models
acf(model_lm4$residuals)
naive_mod = arima(y, order = c(0,1,0), 
      seasonal = list(order = c(0, 1, 0), 
                      period = 12), xreg = xreg)
checkresiduals(naive_mod)
# Modelling ARIMA errors --------------------------------------------------

# Both fail diagnostics 
fit = auto.arima(log_beer_ts, xreg = log_Steel, stepwise = FALSE)
checkresiduals(fit)

fit2 = auto.arima(log_beer_ts, stepwise = FALSE)
checkresiduals(fit2)

# Try regression with the ARIMA model we used in Scenario 3
# Test out some model fits with the new regressors 
fit_reg_1 = arima(
  y,
  order    = c(2, 1, 3),    
  seasonal = c(0, 1, 1)   
  
)

summary(fit_reg_1)
checkresiduals(fit_reg_1)

fit_reg_2 = arima(
  y,
  order = c(2,1,3),
  seasonal = c(0,1,1),
  xreg = xreg
)
checkresiduals(fit_reg_2)


fit_reg_3 = arima(
  y,
  order = c(2,1,3),
  seasonal = c(0,1,1),
  xreg = xreg2
)
checkresiduals(fit_reg_3)

fit_reg_4 = arima(
  y,
  order = c(2,1,3),
  seasonal = c(0,1,1),
  xreg = xreg3
)
checkresiduals(fit_reg_4)

fit_reg_5 = arima(
  y,
  order = c(2,1,3),
  seasonal = c(0,1,1),
  xreg = xreg4
)
checkresiduals(fit_reg_5)

fit_reg_6 = arima(
  y,
  order = c(2,1,3),
  seasonal = c(0,1,1),
  xreg = xreg5
)
checkresiduals(fit_reg_6)


model_reg = list(fit_reg_1, fit_reg_2, fit_reg_3, fit_reg_4, fit_reg_5, fit_reg_6)

# Validation --------------------------------------------------------------
candidate_models = list(
  list(name = "xreg1", desc = "Elec, Gas, Steel", xreg = xreg1),
  list(name = "xreg2", desc = "Gas, Steel", xreg = xreg2),
  list(name = "xreg3", desc = "Elec, Elec^2", xreg = xreg3),
  list(name = "xreg4", desc = "Steel, Elec", xreg = xreg4),
  list(name = "xreg5", desc = "Elec, Gas", xreg = xreg5)
)

# For each of the regressors we will simply take the auto arima as the 
forecast_regressor = function(x, h) {
  forecast(auto.arima(x), h = h)$mean
}

results = data.frame(Model = character(), Description = character(), MSE = numeric(), stringsAsFactors = FALSE)

# Number of rolling validation sets
num_splits = 5

for (candidate in candidate_models) {
  model_name = candidate$name
  desc = candidate$desc
  xreg_full = as.matrix(ts(candidate$xreg, start = start(y), frequency = frequency(y)))
  
  n_total = length(y)
  mse_vals = numeric(num_splits)
  
  for (k in 1:num_splits) {
    end_train_index = n_total - h * (num_splits - k + 1)
    
    train_y = window(y, end = time(y)[end_train_index])
    test_y  = window(y, start = time(y)[end_train_index + 1], end = time(y)[end_train_index + h])
    
    train_xreg = window(xreg_full, end = time(y)[end_train_index])
    fc_xreg = forecast_regressors(train_xreg, h)
    
    fit_model = Arima(train_y, order = c(2,1,3), seasonal = c(0,1,1),
                      xreg = train_xreg, include.constant = TRUE, method = "ML")
    
    fc_y = forecast(fit_model, h = h, xreg = fc_xreg)
    
    acc = accuracy(fc_y, test_y)
    mse_vals[k] = acc["Test set", "RMSE"]^2
    print("here")
  }
  
  avg_mse = mean(mse_vals)
  
  results = rbind(results, data.frame(Model = model_name,
                                      Description = desc,
                                      MSE = avg_mse,
                                      stringsAsFactors = FALSE))
}

results = results %>% arrange(MSE)
print(results)

# Add in the baseline model
mse_vals_baseline = numeric(5)

for (k in 1:num_splits) {
  end_train_index = n_total - h * (num_splits - k + 1)
  
  train_y = window(y, end = time(y)[end_train_index])
  test_y  = window(y, start = time(y)[end_train_index + 1], end = time(y)[end_train_index + h])
  
  fit_baseline = Arima(train_y, order = c(2,1,3), seasonal = c(0,1,1),
                       include.constant = TRUE, method = "ML")
  
  fc_y = forecast(fit_baseline, h = h)
  
  acc = accuracy(fc_y, test_y)
  mse_vals_baseline[k] = acc["Test set", "RMSE"]^2
}

avg_mse_baseline = mean(mse_vals_baseline)

# Add to results table
results = rbind(results, data.frame(Model = "baseline",
                                    Description = "ARIMA(2,1,3)(0,1,1)[12] with no xreg",
                                    MSE = avg_mse_baseline,
                                    stringsAsFactors = FALSE))

# FINAL MODEL --------------------------------------------------------------
# get the predictions for the future 
fc_Gas = forecast_regressor(Gas, h)
fc_Steel = forecast_regressor(Steel, h)

# Fit the final model 
fc_xreg2 = cbind(Gas = fc_Gas, Steel = fc_Steel)
fit_final_4 = Arima(y, order = c(2,1,3), seasonal = c(0,1,1),
                  xreg = xreg2, include.constant = TRUE, method = "ML")

# Seems reasonable
fc_y = forecast(fit_final_4, h = h, xreg = fc_xreg2)
autoplot(fc_y) +
  ggtitle("24-Day-Ahead Forecast using xreg2 (Gas and Steel)") +
  ylab("log(Beer Production)") +
  xlab("Time")

# Transfrom back to original scale
fc_y_exp = fc_y
fc_y_exp$mean = exp(fc_y$mean)
fc_y_exp$lower = exp(fc_y$lower)
fc_y_exp$upper = exp(fc_y$upper)

autoplot(fit2_imp, series = "Beer Production") + 
  autolayer(fc_y_exp, series = "Forecast", PI = TRUE) +
  ggtitle("24-Month-Ahead Forecast of Beer Production") +
  xlab("Time") +
  ylab("Beer Production") +
  scale_color_manual(values = c("Beer Production" = "darkblue", 
                                "Forecast" = "red")) +
  guides(colour = guide_legend(title = "Series"))


