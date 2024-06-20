# code to prepare `df_uc` and `df_uc_source`

library(tidyverse)

set.seed(1234)
n <- 100000
effect_strength <- 2

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
c2 <- rbinom(n, 1, 0.2)
c3 <- rbinom(n, 1, 0.8)
u <- rbinom(n, 1, 0.5)
x_bi <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(0.75) * c2 +
                              log(2.5) * c3 + log(2) * u))
x_cont <- -2 + 1.5 * c1 + 0.75 * c2 + 2.5 * c3 + 2 * u + rnorm(n, 0, 1)
y_bi <- rbinom(n, 1, plogis(-2.5 + log(effect_strength) * x_bi + log(1.5) * c1 -
                              log(2.5) * c2 - log(0.75) * c3 + log(2) * u))
y_cont <- (
  -2.5 + effect_strength * x_cont + 1.5 * c1 - 2.5 * c2 - 0.75 * c3 + 2 * u +
    rnorm(n, 0, 1)
)

df <- data.frame(
  X_bi = x_bi,
  X_cont = round(x_cont, 2),
  Y_bi = y_bi,
  Y_cont = round(y_cont, 2),
  C1 = c1,
  C2 = c2,
  C3 = c3,
  U = u
)

rm(c1, c2, c3, u, x_bi, x_cont, y_bi, y_cont)

# INSPECT MODELS
nobias_model <- glm(Y_bi ~ X_bi + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.00 (1.93, 2.07)

bias_model <- glm(Y_bi ~ X_bi + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 2.22 (2.14, 2.29)

# OBTAIN BIAS PARAMETERS
u_model <- glm(U ~ X_bi + Y_bi + C1 + C2 + C3,
               data = df,
               family = binomial(link = "logit"))

# ADJUST
est <- adjust_uc(
  df,
  "X_bi",
  "Y_bi",
  c("C1", "C2", "C3"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4],
    u_model$coef[5],
    u_model$coef[6]
  )
)

between(est$estimate, 1.9, 2.1)

# CREATE PACKAGE DATA
df_uc_source <- df
head(df_uc_source)
use_data(df_uc_source, overwrite = TRUE)

df_uc <- df %>%
  select(-U) # no access to this in the real world
head(df_uc)
use_data(df_uc, overwrite = TRUE)
