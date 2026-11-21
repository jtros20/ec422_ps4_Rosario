library(tidyverse)
library(fpp3)
library(BVAR)

liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv",
  na = ".",
  col_types = cols(
    observation_date = col_date(format = "%Y-%m-%d"),
    MRTSSM4453USN = col_double()
  )
)

liquor <- liquor_raw |>
  transmute(Month = yearmonth(observation_date), sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month <= yearmonth("2019 Dec"))

has_gaps(liquor)

# Plot liquor sales over time
liquor |>
  autoplot(sales) +
  labs(
    title = "Monthly U.S. Liquor Store Sales",
    x = "Year",
    y = "Sales (Millions of Dollars)"
  ) +
  theme_minimal()

# Season plot of logged liquor sales
liquor |>
  gg_season(log(sales)) +
  labs(
    title = "Seasonal Pattern in U.S. Liquor Store Sales",
    x = "Month",
    y = "Log of Sales"
  ) +
  theme_minimal()

#Sales in U.S. liquor stores are increasing overall, although sales peak in December and trough in January and February. Generally, as sales increase, volatility increases. Using logarithms will reduce this volatility.

liquor_fit <- liquor |>
  model(arima = ARIMA(log(sales)),
        ets   = ETS(log(sales)))

liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()

#ARIMA(1,0,2)(2,1,1)[12] with drift was selected. This means that seasonal differencing has been applied (with order 1) with 12 month periodicity. No regular differencing was required. The ETS model selected is ETS(A,A,A). This means that additive errors, additive trend and additive seasonality have been used.

liquor_fit |> select(arima) |> gg_tsresiduals()

liquor_fit |>
  select(arima) |>
  augment() |>
  features(.innov, ljung_box, lag = 24, dof = 6)

liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))

liquor_cv_fit <- liquor_cv |>
  model(snaive = SNAIVE(sales),
        ets    = ETS(log(sales)),
        arima  = ARIMA(log(sales)))

saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds")

# 228 origins, 348 - 120.

liquor_cv_fc <- liquor_cv_fit |>
  forecast(h = 1)

liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))

# Formatted Table
liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE)) |>
  knitr::kable(digits = 2)

# first Forecast is January 2001, last is December 2019. The ARIMA model produces the most accurate forecasts, with the lowest RMSE (81.43) and MAE (65.65). Both metrics agree that ARIMA performs best, followed by ETS and SNAIVE.
# Diebold-M test
liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = sales),
            by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

forecast::dm.test(liquor_errs$arima, liquor_errs$ets, h = 1, power = 2)

# ARIMA vs. Seasonal Naive
forecast::dm.test(
  liquor_errs$arima,
  liquor_errs$snaive,
  h = 1,
  power = 2
)
# The ARIMA vs. ETS test produced a DM statistic of -2.9062 (p = 0.004042), while ARIMA vs. SNAIVE produced -10.755 (p < 2.2e-16). Both results are statistically significant, so I would choose ARIMA for its better forecast accuracy.


