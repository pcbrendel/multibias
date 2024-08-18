# code to prepare `df_omc` and `df_omc_source`

library(tidyverse)

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
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3, Ystar = ystar)

rm(c1, c2, c3, x, y, ystar)

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

bias_model <- glm(Ystar ~ X + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.40 (1.36, 1.44)

# OBTAIN BIAS PARAMETERS
y_model <- glm(Y ~ X + Ystar + C1 + C2 + C3,
               data = df,
               family = binomial(link = "logit"))
summary(y_model)

# ADJUST

adjust_omc(
  df,
  "X",
  "Ystar",
  c("C1", "C2", "C3"),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4],
    y_model$coef[5],
    y_model$coef[6]
  )
)
# 1.91 (1.83, 1.99)

# CREATE PACKAGE DATA
df_omc_source <- df
head(df_omc_source)
use_data(df_omc_source, overwrite = TRUE)

df_omc <- df %>%
  select(X, Ystar, C1, C2, C3) # only have access to these in real-world
head(df_omc)
use_data(df_omc, overwrite = TRUE)
