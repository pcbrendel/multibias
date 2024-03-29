# code to prepare `df_uc_sel` and `df_uc_sel_source`

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
u <- rbinom(n, 1, .5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(3) * u))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1 + log(3) * u))
s <- rbinom(n, 1, plogis(log(2.5) * x + log(2.5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, U = u, S = s)
s1df <- df[sample(seq_len(n), size = n, replace = TRUE,
                  prob = df$S), ]
rm(c1, u, x, y, s)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.02 (1.95, 2.09)

bias_model <- glm(Y ~ X + C1,
                  family = binomial(link = "logit"),
                  data = s1df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 2.24 (2.17, 2.31)

# OBTAIN BIAS PARAMETERS
u_model <- glm(U ~ X + Y + C1,
               data = df,
               family = binomial(link = "logit"))
s_model <- glm(S ~ X + Y,
               data = df,
               family = binomial(link = "logit"))

# ADJUST
adjust_uc_sel(
  s1df,
  "X",
  "Y",
  "C1",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)
# 2.01 (1.96, 2.07)

# CREATE PACKAGE DATA
df_uc_sel_source <- df
head(df_uc_sel_source)
use_data(df_uc_sel_source)

row.names(s1df) <- NULL
df_uc_sel <- s1df %>%
  select(X, Y, C1) # only have access to these in real-world
head(df_uc_sel)
use_data(df_uc_sel, overwrite = TRUE)
