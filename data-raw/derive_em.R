# code to prepare `df_em` and `df_em_source`

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
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))

df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3, Xstar = xstar)

rm(c1, c2, c3, x, y, xstar)

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

bias_model <- glm(Y ~ Xstar + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.51 (1.45, 1.56)

# OBTAIN BIAS PARAMETERS
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df)
summary(x_model)

# ADJUST
adjust_em(
  df,
  "Xstar",
  "Y",
  c("C1", "C2", "C3"),
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4],
    x_model$coef[5],
    x_model$coef[6]
  )
)
# 1.87 (1.79, 1.94)

# CREATE PACKAGE DATA
df_em_source <- df
head(df_em_source)
use_data(df_em_source, overwrite = TRUE)

df_em <- df %>%
  select(Xstar, Y, C1, C2, C3) # only have access to these in real-world
head(df_em)
use_data(df_em, overwrite = TRUE)
