library(tidyverse)
library(fpp3)
library(BVAR)

# Part 2
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

# Part 3
# Equal Weighted ARIMA and ETS
liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |>
  forecast(h = 1)

liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))

#The combo model performs better than ETS, with a lower RMSE (85.2 vs. 94.4) and MAE (66.8 vs. 72.3). However, ARIMA still performs slightly better than the combo model on both metrics.

cor(liquor_errs$arima, liquor_errs$ets)
#A combination of these forecasts do not perform better. This is likely due to strong correlation between the forecast errors from ARIMA and ETS models (correlation = 0.8783). Further, as ARIMA had significantly better predictions than ETS (p-value = 0.0040), a combination of the two models has resulted in very slightly worse forecasts.

# Part 4
# Read csv
qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |>
  slice(-(1:2))

# Building 8 series as quarterly growth rates in percent
growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date,
            gdp    = growth(GDPC1),     # real GDP
            cons   = growth(PCECC96),   # real consumption
            inv    = growth(GPDIC1),    # real investment
            emp    = growth(PAYEMS),    # payroll employment
            infl   = growth(GDPCTPI),   # GDP price index
            unrate = difference(UNRATE),     # unemployment rate
            ffr    = difference(FEDFUNDS),   # federal funds rate
            gs10   = difference(GS10)) |>    # 10-year Treasury yield
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

y <- qd |> select(-date) |> as.matrix()

nrow(y)
head(y)

# Set up TSCV with first 120 quarters
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]

# OLS VAR
var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)            # columns: y_t, then its lags 1 to p
  M <- ncol(z)                    # number of variables
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))   # constant, latest p lags
  x %*% B
}

fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}

# Compute RMSE
rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))

rmse_by_var(fc_var1)
rmse_by_var(fc_var4)

#VAR(1) produces superior forecasts for 7 out of 8 variables. However, although VAR(4) is a generalisation of VAR(1), VAR(1) produces slightly superior forecasts for inflation. It may be that the smaller model produces better forecasts as it has fewer parameters to estimate. This means less overfitting and lower estimation error.

# Estimate BVAR(4) With minnesota prior
bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))         # constant and 4 lags of each variable
    fc <- rbind(fc, x %*% coef(fit_T))   # forecast of y_{T+1}
  }
  fc
}

set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

# Tighter prior
set.seed(422)
fc_bvar_005 <- bvar_cv(0.05)

var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)

saveRDS(var_results, "results/var_results.rds")

var_results |> knitr::kable(digits = 3)

# 1. BVAR(4) with lambda = 0.3 outperforms VAR(4) for all eight variables because Bayesian shrinkage reduces overfitting and estimation error.

# 2. VAR(1) performs better for five variables, while BVAR(4) with lambda = 0.3 performs better for investment, inflation, and the 10-year Treasury yield.

# 3. As expected when lambda=0.05, increasing shrinkage has a negative impact on performance. The optimal model may have an over-optimistic RMSE due to selection bias as the same set of results was used for model selection.

# Diebold-Mariano Tests

# GDP Growth: VAR(4) vs. BVAR(4), lambda = 0.3
e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

# Inflation: VAR(4) vs. BVAR(4), lambda = 0.3
e_var4_infl <- actual[, "infl"] - fc_var4[, 5]
e_bvar_infl <- actual[, "infl"] - fc_bvar_03[, 5]

forecast::dm.test(e_var4_infl, e_bvar_infl, h = 1, power = 2)

# For the forecasts of GDP, the DM statistic equals 2.1366 and is thus significant (p=0.03468). The DM statistic for forecasts of inflation is 1.4529, and therefore insignificant (p=0.1489). This shows that although the BVAR(4) model performs better in terms of RMSE (see 4.3b), it does not necessarily lead to significant improvements in forecast accuracy.

# Part 5

# transformation
md_raw <- read_csv("data/fredmd_2026-09.csv", show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))

transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x,                      # level
         "2" = d(x),                   # first difference
         "3" = d(d(x)),                # second difference
         "4" = log(x),                 # log
         "5" = 100 * d(log(x)),        # growth rate
         "6" = 100 * d(d(log(x))),     # change in the growth rate
         "7" = d(x / dplyr::lag(x) - 1))
}

# Apply transformation
fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

fm <- fm[, colSums(is.na(fm)) == 0]

dim(fm)

tcode[c("INDPRO", "PAYEMS", "CPIAUCSL", "FEDFUNDS")]
# 1. INDPRO (5): First difference of the log to get monthly change in industrial production.

# 2. PAYEMS (5): First difference of the log to get monthly change in employment.

# 3. CPIAUCSL (6): Second difference of the log to get change in inflation.

# 4. FEDFUNDS (2): First difference of the variable to get monthly change in interest rate.

# 5. Note that FEDFUNDS has been transformed differently. Since interest rates have been given as a percentage we take the first difference. However, in INDPRO we have a positive index and so it makes more sense to take percentage change.

# Estimate factors
Z <- fm |> select(-date) |> as.matrix()

pc <- prcomp(scale(Z))

summary(pc)$importance[, 1:5]

# The first three factors explain 29.45% of the total variation across all 121 series.

# Flip first factor
F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))

tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) +
  geom_line() +
  labs(title = "First Factor from FRED-MD", x = NULL, y = NULL) +
  theme_minimal()

# Eight series most correlated
sort(abs(cor(pc$x[, 1], Z))[1, ], decreasing = TRUE)[1:8]
sort(abs(cor(pc$x[, 2], Z))[1, ], decreasing = TRUE)[1:8]

# 1. The first factor captures overall economic activity, especially employment and industrial production. The large declines during recessions, particularly around 2008–2009, support this interpretation.
# 2. The second factor captures inflation and price movements, since its strongest correlations are with the Consumer Price Index and personal consumption expenditure price measures.

# replicate horse race
ip <- fm$INDPRO
n <- nrow(Z)
init <- which(fm$date == as.Date("1989-12-01"))

e_ar <- NULL
e_fa <- NULL
for (T in init:(n - 1)) {
  F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:3]         # factors using data through T
  d <- data.frame(y_next = ip[2:T], y = ip[1:(T - 1)], F_T[1:(T - 1), ])
  new <- data.frame(y = ip[T], t(F_T[T, ]))
  e_ar <- c(e_ar, ip[T + 1] - predict(lm(y_next ~ y, d), new))
  e_fa <- c(e_fa, ip[T + 1] - predict(lm(y_next ~ ., d), new))
}

c(RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)))

# This matches values from class. The RMSE is 0.617 (AR(1)) and 0.589 (Factor Model). However, using the whole sample to estimate factors gives biased forecasts, as the values in the forecasts would not be known at the time of forecasting, meaning there is look ahead bias.

# Direct forecasting horse race
direct_race <- function(series, h = 12, r = 3) {
  lv <- md |>
    mutate(log_x = log(.data[[series]]),
           g = 100 * (log_x - lag(log_x))) |>
    filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
  yh <- 100 * (lead(lv$log_x, h) - lv$log_x)      # growth over the next h months
  g  <- lv$g                                      # growth over the latest month
  
  e_ar <- NULL
  e_fa <- NULL
  for (T in init:(n - h)) {
    F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:r, drop = FALSE]
    i <- 1:(T - h)
    d <- data.frame(yh = yh[i], g = g[i], F_T[i, , drop = FALSE])
    new <- data.frame(g = g[T], F_T[T, , drop = FALSE])
    e_ar <- c(e_ar, yh[T] - predict(lm(yh ~ g, d), new))
    e_fa <- c(e_fa, yh[T] - predict(lm(yh ~ ., d), new))
  }
  tibble(series, forecasts = length(e_ar),
         RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)),
         MAE_AR = mean(abs(e_ar)), MAE_FA = mean(abs(e_fa)),
         DM_pvalue = forecast::dm.test(e_ar, e_fa, h = h, power = 2,
                                       varestimator = "bartlett")$p.value)
}

direct_race("CPIAUCSL")

# The factors improve 12-month CPI inflation forecasts, reducing RMSE from 2.11 to 1.89, a 10.43% reduction. The improvement is statistically significant, with a DM test p-value of 0.000443.