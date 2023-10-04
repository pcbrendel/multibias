# code to prepare `df_omc` and `df_omc_source`

library(tidyverse)
library(nnet)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, Ystar = ystar)

rm(c1, x, y, ystar)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 1.97 (1.87, 2.07)

bias_model <- glm(Ystar ~ X + C1,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.44 (1.39, 1.49)

# OBTAIN BIAS PARAMETERS
y_model <- glm(Y ~ X + Ystar + C1,
               data = df,
               family = binomial(link = "logit"))
summary(y_model)

# ADJUST

adjust_omc(
  df,
  "X",
  "Ystar",
  "C1",
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4]
  )
)
# 1.98 (1.88, 2.08)

# CREATE PACKAGE DATA
df_omc_source <- df
head(df_omc_source)
use_data(df_omc_source)

df_omc <- df %>%
  select(X, Ystar, C1) # only have access to these in real-world
head(df_omc)
use_data(df_omc)
