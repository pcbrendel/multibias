library(devtools)
library(tidyverse)
library(roxygen2)

# library(urlchecker)
# library(gitcreds)
# library(available)

# # registerS3method("print", "data_observed", print.data_observed)
# # registerS3method("print", "data_validation", print.data_validation)

load_all()
build()

# document
roxygenize()
document()

# testing
usethis::use_testthat()
use_test("adjust_em_sel") # creates test
test_file("tests/testthat/test-adjust_uc.R") # single test
test() # tests all
devtools::run_examples(".")

# check
check()
check(remote = TRUE, manual = TRUE)
devtools::check_win_devel()

# submit to cran
usethis::use_release_issue()
usethis::use_version("minor")
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
usethis::use_roxygen_md()
usethis::use_lifecycle()

# vignette
usethis::use_vignette("my-vignette")
devtools::build_rmd("vignettes/multibias.Rmd")

# git PAT
create_github_token()
gh_token_help()
gitcreds::gitcreds_set()
git_sitrep()


# README example 2

head(df_uc_em_sel)

biased_model <- glm(Y ~ Xstar + C1 + C2 + C3,
  data = df_uc_emc_sel,
  family = binomial(link = "logit")
)
biased_or <- round(exp(coef(biased_model)[2]), 2)
print(paste0("Biased Odds Ratio: ", biased_or))

u_model <- glm(U ~ X + Y,
  data = df_uc_emc_sel_source,
  family = binomial(link = "logit")
)
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
  data = df_uc_emc_sel_source,
  family = binomial(link = "logit")
)
s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
  data = df_uc_emc_sel_source,
  family = binomial(link = "logit")
)

library(doParallel)

no_cores <- detectCores() - 1
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores)

set.seed(1234)
nreps <- 1000
est <- vector(length = nreps)

or <- foreach(
  i = 1:nreps, .combine = c,
  .packages = c("dplyr", "multibias")
) %dopar% {
  df_sample <- df_uc_em_sel[sample(
    seq_len(nrow(df_uc_emc_sel)), nrow(df_uc_emc_sel),
    replace = TRUE
  ), ]

  est[i] <- adjust_uc_em_sel(
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
