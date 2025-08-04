# Forecasting Challenge – University of Waterloo (2025)

This repository includes my submission to a graduate forecasting competition held at the University of Waterloo. The task involved developing forecasting models for five different real-world time series. Each scenario posed unique challenges in modeling, data exploration, and evaluation.

Scenario 5 (Pollution Forecasting) was selected as a winning submission. I used dynamic harmonic regression to capture complex seasonal patterns, supported by diagnostics and visualizations.

---

## Files

- `Model_Fitting_Scen1.R` — Scenario 1
- `Model_Fitting_Scen2.R` — Scenario 2
- `Model_Fitting_Scen3.R` — Scenario 3
- `Model_Fitting_Scen4.R` — Scenario 4
- `Model_Fitting_Scen5.R` — Scenario 5 (Pollution)
- `report.pdf` — Full write-up summarizing methods and results

---

## Methods & Tools

Most of the work was done in R using packages like `fpp3`, `forecast`, and `ggplot2`. Techniques include:

- SARIMA, SARIMAX
- ARMA-GARCH
- Kalman filtering
- Dynamic harmonic regression
- Cross-validation for model tuning

Forecast accuracy was assessed using RMSE, MAE, and visual checks for fit and residuals.

---

## Recognition

- Scenario 5 selected as a winning forecast
- Recipient of the Time Series and Forecasting Award from the Business and Industrial Statistics Research Group (BISRG)

---

Feel free to explore the code or report for details on the models and results.
