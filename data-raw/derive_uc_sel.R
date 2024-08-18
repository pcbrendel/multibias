# code to prepare `df_uc_sel` and `df_uc_sel_source`

library(tidyverse)

set.seed(1234)
n <- 100000
effect_strength <- 2

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
c2 <- rbinom(n, 1, 0.2)
c3 <- rbinom(n, 1, 0.8)
u <- rbinom(n, 1, .5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(0.75) * c2 +
                           log(2.5) * c3 + log(2) * u))
y <- rbinom(n, 1, plogis(-2.5 + log(effect_strength) * x + log(1.5) * c1 -
                           log(2.5) * c2 - log(0.75) * c3 + log(2) * u))
s <- rbinom(n, 1, plogis(log(2.5) * x + log(2.5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3, U = u, S = s)
s1df <- df[sample(seq_len(n), size = n, replace = TRUE,
                  prob = df$S), ]
rm(c1, c2, c3, u, x, y, s)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 1.99 (1.93, 2.07)

bias_model <- glm(Y ~ X + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = s1df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.87 (1.81, 1.93)

# OBTAIN BIAS PARAMETERS
u_model <- glm(U ~ X + Y + C1 + C2 + C3,
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
  c("C1", "C2", "C3"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4],
    u_model$coef[5],
    u_model$coef[6]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)
# 2.00 (1.95, 2.06)

# CREATE PACKAGE DATA
df_uc_sel_source <- df
head(df_uc_sel_source)
use_data(df_uc_sel_source, overwrite = TRUE)

row.names(s1df) <- NULL
df_uc_sel <- s1df %>%
  select(X, Y, C1, C2, C3) # only have access to these in real-world
head(df_uc_sel)
use_data(df_uc_sel, overwrite = TRUE)
