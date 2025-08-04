library(astsa)
library(fpp3) 
library(tidyverse)
library(forecast)
library(imputeTS)

prod_targ = read.csv("Data/ProdData/prod_target.txt")
prod_target_ts = ts(prod_targ[,2])
prod_targ <- prod_targ %>% 
  select(V1, V2) %>% 
  rename(time = V1, Beer = V2) %>% 
  mutate(time = yearmonth(time)) %>% # This is a random time.
  as_tsibble(index = time)
            
# Look at this plot. What an amazing package     
prod_targ %>% gg_tsdisplay() + labs(title = "Monthly Beer Production")

prod_targ %>% gg_tsdisplay(difference(Beer, 6) , 
                           plot_type = c("partial"), lag_max = 60)
prod_targ %>% gg_tsdisplay(difference(Beer, 12) , 
                           plot_type = c("partial"), lag_max = 60)


beer_ts <- ts(prod_targ$Beer, frequency = 12)



# Model Fitting Algo for State Space --------------------------------------

results <- list()
model_num <- 1

# Grid search parameters
p_vals <- 0:3
d_vals <- 0:1
q_vals <- 0:3
P_vals <- 0:3
Q_vals <- 0:3
D <- 1
seasonal_periods <- c(6, 12)

# Grid search
for (s in seasonal_periods) {
  for (d in d_vals) {
    for (p in p_vals) {
      for (q in q_vals) {
        for (P in P_vals) {
          for (Q in Q_vals) {
            try({
              fit <- Arima(beer_ts,
                           order = c(p, d, q),
                           seasonal = list(order = c(P, D, Q), period = s),
                           include.constant = FALSE)
              
              lb_pval <- Box.test(residuals(fit), lag = 24,
                                  fitdf = p + q + P + Q, type = "Ljung")$p.value
              
              results[[model_num]] <- list(
                model = fit,
                order = c(p, d, q),
                seasonal = c(P, D, Q),
                seasonal_period = s,
                AIC = AIC(fit),
                LjungBox_p = lb_pval
              )
              model_num <- model_num + 1
            }, silent = TRUE)
          }
        }
      }
    }
  }
}


results_df <- bind_rows(lapply(results, function(x) {
  tibble(
    order = paste0("(", paste(x$order, collapse=","), ")"),
    seasonal = paste0("(", paste(x$seasonal, collapse=","), ")[", x$seasonal_period, "]"),
    AIC = x$AIC,
    LjungBox_p = x$LjungBox_p
  )
}))


# Fits from the model fitting algorithm -----------------------------------

fit1 <- Arima(beer_ts, order = c(2, 1, 3),
              seasonal = list(order = c(1, 1, 1), period = 6),
              include.constant = FALSE)
checkresiduals(fit1)

# These other two fits had a reasoanable fit 
fit2 <- Arima(beer_ts, order = c(2, 1, 3),
              seasonal = list(order = c(0, 1, 1), period = 12),
              include.constant = FALSE)
checkresiduals(fit2)
fit3 <- Arima(beer_ts, order = c(3, 1, 3),
              seasonal = list(order = c(0, 1, 1), period = 12),
              include.constant = FALSE)
checkresiduals(fit3)

# Kalman Smoothing --------------------------------------------------------

# Basic attempt
fit1_imp = na_kalman(beer_ts)
tsplot(fit1_imp)

# Fit with my own state space model from above 
mod_ss = arima(beer_ts, order = c(2,1,3), 
               seasonal = list(order = c(0, 1, 1), 
                               period = 12))


fit2_imp = na_kalman(beer_ts, model = mod_ss$model, smooth = TRUE)
ggplot_na_imputations(beer_ts, fit2_imp)


# Try other imputation algorithms 
fit3_imp <- na_interpolation(prod_target_ts, option = "stine")
ggplot_na_imputations(beer_ts,fit3_imp) + labs(title = "Imputation from Interpoltaion")

fit4_imp <- na_kalman(prod_target_ts, model = "auto.arima")
ggplot_na_imputations(beer_ts,fit4_imp) + labs(title = "Imputation from Kalman Auto Arima")

fit5_imp <- na_ma(prod_target_ts)
ggplot_na_imputations(beer_ts,fit5_imp) + labs(title = "Imputation from Moving Average")


# Confidence intervals ----------------------------------------------------

resid_ss = as.vector(fit2$residuals)
sd_imp = sd(resid_ss, na.rm = TRUE)

# Get the value and index of imputed values 
miss_idx = which(is.na(as.vector(beer_ts)))

imp_val = fit2_imp[miss_idx]
CI_up = imp_val + 2*sd_imp
CI_lw = imp_val - 2*sd_imp



plot(as.vector(beer_ts), type = "l", ylab = "Beer Production", xlab = "Time")
points(miss_idx, imp_val, col = "red", type = "o")
points(miss_idx, CI_up, col = "blue", type = "l")
points(miss_idx, CI_lw, col = "blue", type = "l")
legend("topleft", legend = c("Imputed Value", "95% Confidence Interval"), 
       fill= c("red", "blue"))


# Visual  -----------------------------------------------------------------
library(ggplot2)

# Set up vectors
beer_vec <- as.vector(beer_ts)
time_idx <- seq_along(beer_vec)

# Make a df for ggplot 
df <- data.frame(
  Time = time_idx,
  Beer = beer_vec
)

# Get the imputations
imp_df <- data.frame(
  Time = miss_idx,
  Imputed = imp_val,
  CI_Lower = CI_lw,
  CI_Upper = CI_up
)

ggplot(df, aes(x = Time, y = Beer)) +
  geom_line(color = "black") +
  geom_ribbon(data = imp_df, aes(x = Time, ymin = CI_Lower, ymax = CI_Upper), 
              fill = "blue", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = imp_df, aes(x = Time, y = Imputed), color = "red", size = 0.5) +
  geom_point(data = imp_df, aes(x = Time, y = Imputed), color = "red", size = 1.2) +
  labs(title = "Beer Production with Imputed Values",
       x = "Time", y = "Beer Production") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) 
