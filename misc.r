library(devtools)
library(tidyverse)
library(available)
library(roxygen2)
library(urlchecker)

install.packages("gitcreds")

# submit to cran
usethis::use_release_issue()
usethis::use_version("major")
devtools::submit_cran()

# other
available("multibias")
create_package("C:/Users/brend/Desktop/folder/multibias")
use_data_raw()
use_mit_license()
use_github_action()
use_lifecycle() # when functions become experimental/superseded/deprecated
use_news_md()
use_cran_comments()
urlchecker::url_check()
usethis::use_package("rlang", min_version = TRUE)

# git PAT
create_github_token()
gh_token_help()
gitcreds::gitcreds_set()
git_sitrep()

# document
roxygenize()
document()

# testing
usethis::use_testthat()
use_test("adjust_emc_sel")
use_test("adjust_emc")
use_test("adjust_multinom_uc_emc_sel")
use_test("adjust_multinom_uc_emc")
use_test("adjust_sel")
use_test("adjust_uc_emc_sel")
use_test("adjust_uc_emc")
use_test("adjust_uc_sel")
use_test("adjust_uc")

test()
check()
devtools::check(remote = TRUE, manual = TRUE)
devtools::check_win_devel()

# informal run
load_all()

adjust_sel(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "HPT",
  s_model_coefs = c(qlogis(0.25), log(0.75), log(0.75))
)

# README example 1

head(evans)

biased_model <- glm(CHD ~ SMK + HPT,
                    data = evans,
                    family = binomial(link = "logit"))
or <- round(exp(coef(biased_model)[2]), 2)
or_ci_low <- round(exp(coef(biased_model)[2] -
                         1.96 * summary(biased_model)$coef[2, 2]), 2)
or_ci_high <- round(exp(coef(biased_model)[2] +
                          1.96 * summary(biased_model)$coef[2, 2]), 2)

print(paste0("Biased Odds Ratio: ", or))
print(paste0("95% CI: (", or_ci_low, ", ", or_ci_high, ")"))

mean(evans[evans$SMK == 1, "AGE"])
mean(evans[evans$SMK == 0, "AGE"])

mean(evans[evans$CHD == 1, "AGE"])
mean(evans[evans$CHD == 0, "AGE"])

ggplot(evans, aes(y = factor(SMK), x = AGE, fill = factor(CHD))) + geom_violin()

evans$AGE_bin <- if_else(evans$AGE >= 60, 1, 0)

u_model <- glm(AGE_bin ~ SMK + CHD + HPT,
               data = evans,
               family = binomial(link = "logit"))
summary(u_model)

# int: 25% probability someone is >60 when they are a non-smoker and have no
#      CHD or HPT
# SMK: OR=0.5
# CHD: OR=2.5
# HPT: OR=2 
u_0 <- qlogis(0.25)
u_x <- log(0.5)
u_y <- log(2.5)
u_c <- log(2)

adjust_uc(
  data = evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "HPT",
  u_model_coefs = c(u_0, u_x, u_y, u_c)
)
# 2.27 (1.25, 4.10)

full_model <- glm(CHD ~ SMK + HPT + AGE,
                    data = evans,
                    family = binomial(link = "logit"))
or <- round(exp(coef(full_model)[2]), 2)
or_ci_low <- round(exp(coef(biased_model)[2] -
                         1.96 * summary(full_model)$coef[2, 2]), 2)
or_ci_high <- round(exp(coef(biased_model)[2] +
                          1.96 * summary(full_model)$coef[2, 2]), 2)

print(paste0("Odds Ratio: ", or))
print(paste0("95% CI: (", or_ci_low, ", ", or_ci_high, ")"))


# README example 2

head(df_uc_emc_sel)

biased_model <- glm(Y ~ Xstar + C1 + C2 + C3,
                    data = df_uc_emc_sel,
                    family = binomial(link = "logit"))
biased_or <- round(exp(coef(biased_model)[2]), 2)
print(paste0("Biased Odds Ratio: ", biased_or))

u_model <- glm(U ~ X + Y,
               data = df_uc_emc_sel_source,
               family = binomial(link = "logit"))
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
               data = df_uc_emc_sel_source,
               family = binomial(link = "logit"))
s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
               data = df_uc_emc_sel_source,
               family = binomial(link = "logit"))

library(doParallel)

no_cores <- detectCores() - 1
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores)

set.seed(1234)
nreps <- 1000
est <- vector(length = nreps)

or <- foreach(i = 1:nreps, .combine = c,
              .packages = c("dplyr", "multibias")) %dopar% {

  df_sample <- df_uc_emc_sel[sample(seq_len(nrow(df_uc_emc_sel)),
                                    nrow(df_uc_emc_sel),
                                    replace = TRUE), ]

  est[i] <- adjust_uc_emc_sel(
    df_sample,
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    u_model_coefs = c(
      rnorm(1, mean = u_model$coef[1], sd = summary(u_model)$coef[1, 2]),
      rnorm(1, mean = u_model$coef[2], sd = summary(u_model)$coef[2, 2]),
      rnorm(1, mean = u_model$coef[3], sd = summary(u_model)$coef[3, 2])
    ),
    x_model_coefs = c(
      rnorm(1, mean = x_model$coef[1], sd = summary(x_model)$coef[1, 2]),
      rnorm(1, mean = x_model$coef[2], sd = summary(x_model)$coef[2, 2]),
      rnorm(1, mean = x_model$coef[3], sd = summary(x_model)$coef[3, 2]),
      rnorm(1, mean = x_model$coef[4], sd = summary(x_model)$coef[4, 2]),
      rnorm(1, mean = x_model$coef[5], sd = summary(x_model)$coef[5, 2]),
      rnorm(1, mean = x_model$coef[6], sd = summary(x_model)$coef[6, 2])
    ),
    s_model_coefs = c(
      rnorm(1, mean = s_model$coef[1], sd = summary(s_model)$coef[1, 2]),
      rnorm(1, mean = s_model$coef[2], sd = summary(s_model)$coef[2, 2]),
      rnorm(1, mean = s_model$coef[3], sd = summary(s_model)$coef[3, 2]),
      rnorm(1, mean = s_model$coef[4], sd = summary(s_model)$coef[4, 2]),
      rnorm(1, mean = s_model$coef[5], sd = summary(s_model)$coef[5, 2]),
      rnorm(1, mean = s_model$coef[6], sd = summary(s_model)$coef[6, 2])
    )
  )[[1]]
}

round(median(or), 2)
round(quantile(or, c(.025, .975)), 2)
