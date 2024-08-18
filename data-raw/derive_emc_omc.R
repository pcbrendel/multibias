# code to prepare `df_emc_omc` and `df_emc_omc_source`

library(tidyverse)
library(nnet)

set.seed(1234)
n <- 100000
effect_strength <- 2

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
c2 <- rbinom(n, 1, 0.2)
c3 <- rbinom(n, 1, 0.8)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(0.75) * c2 +
                           log(2.5) * c3))
y <- rbinom(n, 1, plogis(-2.5 + log(effect_strength) * x + log(1.5) * c1 -
                           log(2.5) * c2 - log(0.75) * c3))
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3,
                 Xstar = xstar, Ystar = ystar)

rm(c1, c2, c3, x, y, xstar, ystar)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1 + C2 + C3,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 1.94 (1.86, 2.02)

bias_model <- glm(Ystar ~ Xstar + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.16 (1.13, 1.20)

# OBTAIN BIAS PARAMETERS
# x_model <- glm(X ~ C1 + C2 + C3,
#                data = df,
#                family = binomial(link = "logit"))
# y_model <- glm(Y ~ C1 + C2 + C3,
#                data = df,
#                family = binomial(link = "logit"))

xy_model <- multinom(
  paste0(X, Y) ~ Xstar + Ystar + C1 + C2 + C3,
  data = df
)
summary(xy_model)

# ADJUST
# adjust_emc_omc(
#   df,
#   "Xstar",
#   "Ystar",
#   c("C1", "C2", "C3"),
#   x_model_coefs = c(
#     x_model$coef[1],
#     x_model$coef[2]
#   ),
#   y_model_coefs = c(
#     y_model$coef[1],
#     y_model$coef[2]
#   )
# )
# (, )

adjust_multinom_emc_omc(
  df,
  "Xstar",
  "Ystar",
  c("C1", "C2", "C3"),
  x1y0_model_coefs = c(
    summary(xy_model)$coefficients[2, 1],
    summary(xy_model)$coefficients[2, 2],
    summary(xy_model)$coefficients[2, 3],
    summary(xy_model)$coefficients[2, 4],
    summary(xy_model)$coefficients[2, 5],
    summary(xy_model)$coefficients[2, 6]
  ),
  x0y1_model_coefs = c(
    summary(xy_model)$coefficients[1, 1],
    summary(xy_model)$coefficients[1, 2],
    summary(xy_model)$coefficients[1, 3],
    summary(xy_model)$coefficients[1, 4],
    summary(xy_model)$coefficients[1, 5],
    summary(xy_model)$coefficients[1, 6]
  ),
  x1y1_model_coefs = c(
    summary(xy_model)$coefficients[3, 1],
    summary(xy_model)$coefficients[3, 2],
    summary(xy_model)$coefficients[3, 3],
    summary(xy_model)$coefficients[3, 4],
    summary(xy_model)$coefficients[3, 5],
    summary(xy_model)$coefficients[3, 6]
  )
)
# 1.94 (1.86, 2.02)

# CREATE PACKAGE DATA
df_emc_omc_source <- df
head(df_emc_omc_source)
use_data(df_emc_omc_source, overwrite = TRUE)

df_emc_omc <- df %>%
  select(Xstar, Ystar, C1, C2, C3) # only have access to these in real-world
head(df_emc_omc)
use_data(df_emc_omc, overwrite = TRUE)
