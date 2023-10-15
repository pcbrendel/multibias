# code to prepare `df_uc_omc` and `df_uc_omc_source`

library(tidyverse)
library(nnet)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
u <- rbinom(n, 1, .5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(2) * u))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1 + log(2) * u))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, U = u, Ystar = ystar)

rm(c1, u, x, y, ystar)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.03 (1.95, 2.11)

bias_model <- glm(Ystar ~ X + C1,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.48 (1.43, 1.53)

# OBTAIN BIAS PARAMETERS
u_model <- glm(U ~ X + Y,
               data = df,
               family = binomial(link = "logit"))
y_model <- glm(Y ~ X + Ystar + C1,
               data = df,
               family = binomial(link = "logit"))

uy_model <- multinom(
  paste0(U, Y) ~ X + Ystar + C1,
  data = df
)
summary(uy_model)

# ADJUST
adjust_uc_omc(
  df,
  "X",
  "Ystar",
  "C1",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4]
  )
)
# 2.03 (1.96, 2.11)

adjust_multinom_uc_omc(
  df,
  "X",
  "Ystar",
  "C1",
  u1y0_model_coefs = c(
    summary(uy_model)$coefficients[2, 1],
    summary(uy_model)$coefficients[2, 2],
    summary(uy_model)$coefficients[2, 3],
    summary(uy_model)$coefficients[2, 4]
  ),
  u0y1_model_coefs = c(
    summary(uy_model)$coefficients[1, 1],
    summary(uy_model)$coefficients[1, 2],
    summary(uy_model)$coefficients[1, 3],
    summary(uy_model)$coefficients[1, 4]
  ),
  u1y1_model_coefs = c(
    summary(uy_model)$coefficients[3, 1],
    summary(uy_model)$coefficients[3, 2],
    summary(uy_model)$coefficients[3, 3],
    summary(uy_model)$coefficients[3, 4]
  )
)
# 2.03 (1.95, 2.11)

# CREATE PACKAGE DATA
df_uc_omc_source <- df
head(df_uc_omc_source)
use_data(df_uc_omc_source)

df_uc_omc <- df %>%
  select(X, Ystar, C1) # only have access to these in real-world
head(df_uc_omc)
use_data(df_uc_omc)
