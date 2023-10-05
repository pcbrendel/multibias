# code to prepare `df_uc` and `df_uc_source`

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
u <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(2) * u))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1 + log(2) * u))

df <- data.frame(X = x, Y = y, C1 = c1, U = u)

rm(c1, u, x, y)

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

bias_model <- glm(Y ~ X + C1,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 2.25 (2.16, 2.34)

# OBTAIN BIAS PARAMETERS
u_model <- glm(U ~ X + Y + C1,
               data = df,
               family = binomial(link = "logit"))

# ADJUST
adjust_uc(
  df,
  "X",
  "Y",
  "C1",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4]
  )
)
# 2.02 (1.95, 2.10)

# CREATE PACKAGE DATA
df_uc_source <- df
head(df_uc_source)
use_data(df_uc_source)

df_uc <- df %>%
  select(X, Y, C1) # only have access to these in real-world
head(df_uc)
use_data(df_uc)
