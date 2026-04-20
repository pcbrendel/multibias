library(devtools)
library(tidyverse)
library(roxygen2)

library(usethis)
library(pkgdown)
# library(urlchecker)
# library(gitcreds)
# library(available)

# pkgdown
# usethis::use_pkgdown_github_pages()
build_readme() # sync the .md with the .Rmd
build_site() # local version

# pandoc
rmarkdown::pandoc_version()
rmarkdown::pandoc_available()
rmarkdown::find_pandoc()
usethis::edit_r_environ()
Sys.getenv("RSTUDIO_PANDOC")
Sys.setenv(RSTUDIO_PANDOC = "/Users/pbrendel/.local/share/pandoc")

load_all()
build()

# document
roxygenize()
document()

# testing
usethis::use_testthat()
use_test("adjust_em_sel") # creates test
test_active_file("tests/testthat/test-adjust_uc_om_sel.R") # single test
test() # tests all
devtools::run_examples(".")

# check
check()
check(vignettes = FALSE)
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

# ----

# df_observed <- data_observed(
#   data = df_em,
#   bias = "em",
#   exposure = "Xstar",
#   outcome = "Y",
#   confounders = "C1"
# )

# bp1 <- bias_params(coef_list = list(x = c(-2.10, 1.62, 0.63, 0.35)))
# bp2 <- bias_params(coef_list = list(x = c(-2.10 * 2, 1.62 * 2, 0.63 * 2, 0.35 * 2)))

# result1 <- multibias_adjust(
#   data_observed = df_observed,
#   bias_params = bp1
# )
# result2 <- multibias_adjust(
#   data_observed = df_observed,
#   bias_params = bp2
# )

# multibias_plot(
#   data_observed = df_observed,
#   multibias_result_list = list(
#     "Adjusted with bias params" = result1,
#     "Adjusted with bias params doubled" = result2
#   )
# )
