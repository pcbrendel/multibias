# code to prepare `df_uc_emc` and `df_uc_emc_source`

library(tidyverse)
library(nnet)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
c2 <- rbinom(n, 1, 0.2)
c3 <- rbinom(n, 1, 0.8)
u <- rbinom(n, 1, .5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(0.75) * c2 + log(2.5) * c3 +
                           log(2) * u))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1 - log(2.5) * c2 -
                           log(0.75) * c3))
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))

df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3, U = u, Xstar = xstar)

rm(c1, c2, c3, u, x, y, xstar)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.02 (1.94, 2.09)

bias_model <- glm(Y ~ Xstar + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.53 (1.47, 1.59)

# OBTAIN BIAS PARAMETERS
u_model <- glm(U ~ X + Y,
               data = df,
               family = binomial(link = "logit"))
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
               data = df,
               family = binomial(link = "logit"))

xu_model <- multinom(
  paste0(X, U) ~ Xstar + Y + C1 + C2 + C3,
  data = df
)
summary(xu_model)

# ADJUST
adjust_uc_emc(
  df,
  "Xstar",
  "Y",
  c("C1", "C2", "C3"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4],
    x_model$coef[5],
    x_model$coef[6]
  )
)
# 2.04 (1.96, 2.12)

adjust_multinom_uc_emc(
  df,
  "Xstar",
  "Y",
  c("C1", "C2", "C3"),
  x1u0_model_coefs = c(
    summary(xu_model)$coefficients[2, 1],
    summary(xu_model)$coefficients[2, 2],
    summary(xu_model)$coefficients[2, 3],
    summary(xu_model)$coefficients[2, 4],
    summary(xu_model)$coefficients[2, 5],
    summary(xu_model)$coefficients[2, 6]
  ),
  x0u1_model_coefs = c(
    summary(xu_model)$coefficients[1, 1],
    summary(xu_model)$coefficients[1, 2],
    summary(xu_model)$coefficients[1, 3],
    summary(xu_model)$coefficients[1, 4],
    summary(xu_model)$coefficients[1, 5],
    summary(xu_model)$coefficients[1, 6]
  ),
  x1u1_model_coefs = c(
    summary(xu_model)$coefficients[3, 1],
    summary(xu_model)$coefficients[3, 2],
    summary(xu_model)$coefficients[3, 3],
    summary(xu_model)$coefficients[3, 4],
    summary(xu_model)$coefficients[3, 5],
    summary(xu_model)$coefficients[3, 6]
  )
)
# 2.01 (1.94, 2.09)

# CREATE PACKAGE DATA
df_uc_emc_source <- df
head(df_uc_emc_source)
use_data(df_uc_emc_source, overwrite = TRUE)

df_uc_emc <- df %>%
  select(Xstar, Y, C1, C2, C3) # only have access to these in real-world
head(df_uc_emc)
use_data(df_uc_emc, overwrite = TRUE)