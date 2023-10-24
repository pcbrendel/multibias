# code to prepare `df_omc_sel` and `df_omc_sel_source`

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))
s <- rbinom(n, 1, plogis(log(2.5) * x + log(2.5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, Ystar = ystar, S = s)
s1df <- df[sample(seq_len(n), size = n, replace = TRUE,
                  prob = df$S), ]
rm(c1, x, y, ystar, s)

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
                  data = s1df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.39 (1.34, 1.43)

# OBTAIN BIAS PARAMETERS
y_model <- glm(Y ~ X + Ystar + C1,
               data = df,
               family = binomial(link = "logit"))
s_model <- glm(S ~ X + Ystar + C1,
               data = df,
               family = binomial(link = "logit"))

# ADJUST
adjust_omc_sel(
  s1df,
  "X",
  "Ystar",
  "C1",
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4]
  )
)
# 1.96 (1.89, 2.03)

# CREATE PACKAGE DATA
df_omc_sel_source <- df
head(df_omc_sel_source)
use_data(df_omc_sel_source)

row.names(s1df) <- NULL
df_omc_sel <- s1df %>%
  select(X, Ystar, C1) # only have access to these in real-world
head(df_omc_sel)
use_data(df_omc_sel)