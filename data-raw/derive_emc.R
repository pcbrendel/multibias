# code to prepare `df_emc` and `df_emc_source`

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1))
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))

df <- data.frame(X = x, Y = y, C1 = c1, Xstar = xstar)

rm(c1, x, y, xstar)

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

bias_model <- glm(Y ~ Xstar + C1,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.49 (1.43, 1.55)

# OBTAIN BIAS PARAMETERS
x_model <- glm(X ~ Xstar + Y + C1,
               family = binomial(link = "logit"),
               data = df)
summary(x_model)

# ADJUST
adjust_emc(
  df,
  "Xstar",
  "Y",
  "C1",
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4]
  )
)
# 2.00 (1.91, 2.11)

# CREATE PACKAGE DATA
df_emc_source <- df
head(df_emc_source)
use_data(df_emc_source)

df_emc <- df %>%
  select(Xstar, Y, C1) # only have access to these in real-world
head(df_emc)
use_data(df_emc)
