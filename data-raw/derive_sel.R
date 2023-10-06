# code to prepare `df_sel` and `df_sel_source`

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1))
s <- rbinom(n, 1, plogis(log(2.5) * x + log(2.5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, S = s)
s1df <- df[sample(seq_len(n), size = n, replace = TRUE, prob = df$S), ]

rm(c1, x, y, s)

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

bias_model <- glm(Y ~ X + C1,
                  family = binomial(link = "logit"),
                  data = s1df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.63 (1.57, 1.70)

# OBTAIN BIAS PARAMETERS
s_model <- glm(S ~ X + Y,
               data = df,
               family = binomial(link = "logit"))

# ADJUST
adjust_sel(
  s1df,
  "X",
  "Y",
  "C1",
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)
# 1.93 (1.86, 2.00)

# CREATE PACKAGE DATA
df_sel_source <- df
head(df_sel_source)
use_data(df_sel_source)

row.names(s1df) <- NULL
df_sel <- s1df %>%
  select(X, Y, C1) # only have access to these in real-world
head(df_sel)
use_data(df_sel, overwrite = TRUE)