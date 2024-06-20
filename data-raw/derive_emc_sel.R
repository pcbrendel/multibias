# code to prepare `df_emc_sel` and `df_emc_sel_source`

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c1 <- rbinom(n, 1, 0.5)
c2 <- rbinom(n, 1, 0.2)
c3 <- rbinom(n, 1, 0.8)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c1 + log(0.75) * c2 +
                           log(2.5) * c3))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c1 -
                           log(2.5) * c2 - log(0.75) * c3))
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))
s <- rbinom(n, 1, plogis(log(2.5) * x + log(2.5) * y))

df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3, Xstar = xstar, S = s)
s1df <- df[sample(seq_len(n), size = n, replace = TRUE,
                  prob = df$S), ]
rm(c1, c2, c3, x, y, xstar, s)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C1 + C2 + C3,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 1.94 (1.86, 2.01)

bias_model <- glm(Y ~ Xstar + C1 + C2 + C3,
                  family = binomial(link = "logit"),
                  data = s1df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.43 (1.38, 1.48)

# OBTAIN BIAS PARAMETERS
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
               data = df,
               family = binomial(link = "logit"))
s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
               data = df,
               family = binomial(link = "logit"))

# ADJUST
adjust_emc_sel(
  s1df,
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
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4],
    s_model$coef[5],
    s_model$coef[6]
  )
)
# 1.94 (1.88, 2.00)

# CREATE PACKAGE DATA
df_emc_sel_source <- df
head(df_emc_sel_source)
use_data(df_emc_sel_source, overwrite = TRUE)

row.names(s1df) <- NULL
df_emc_sel <- s1df %>%
  select(Xstar, Y, C1, C2, C3) # only have access to these in real-world
head(df_emc_sel)
use_data(df_emc_sel, overwrite = TRUE)