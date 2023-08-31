# code to prepare `df_uc_emis_sel`
# this data corresponds to the simA data in IJE paper

set.seed(1234)
n <- 100000
C <- rbinom(n, 1, .5)
U <- rbinom(n, 1, .5)

X  <- rbinom(n, 1, expit(-2 + log(1.5) * C + log(2) * U))
# C=1, U=1 -> p = .289
# C=0, U=0 -> p = .119

Y  <- rbinom(n, 1, expit(-2.5 + log(2) * X + log(1.5) * C + log(2) * U))
# X=1, C=1, U=1 -> p = .330
# X=0, C=0, U=0 -> p = .0759

Xstar <- rbinom(n, 1, expit(-1 + log(5) * X + log(1.25) * Y))
# X=1, Y=1 -> p = .697
# X=0, Y=0 -> p = .269

S <- rbinom(n, 1, expit(log(2) * X + log(2) * Y))
# X=1, Y=1 -> p = .800
# X=0, Y=0 -> p = .500

df <- data.frame(C, U, X, Xstar, Y, S)
s1df <- df[sample(1:nrow(df), size = n, replace = TRUE, prob = df$S),]
rm(C, U, X, Y, Xstar, S)

# Inspect biased model vs bias-free model

nobias_model <- glm(Y ~ X + C + U, family = binomial(link = "logit"), data = df)
exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] + summary(nobias_model)$coef[2, 2] * qnorm(.025)), 
  exp(summary(nobias_model)$coef[2, 1] + summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.03 (1.95, 2.11)

bias_model <- glm(Y ~ Xstar + C, family = binomial(link = "logit"), data = s1df)
exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] + summary(bias_model)$coef[2, 2] * qnorm(.025)), 
  exp(summary(bias_model)$coef[2, 1] + summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.46 (1.41, 1.50)

#

usethis::use_data(DATASET, overwrite = TRUE)
