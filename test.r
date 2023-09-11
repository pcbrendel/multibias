library(devtools)
library(tidyverse)
library(available)
library(roxygen2)

install.packages("available")

use_data_raw()
use_mit_license()

# check package with 'available'
available("multibias")

# document
roxygenize()
document()

# check and create the package
check()
create_package("C:/Users/brend/Desktop/folder/multibias")

# release to CRAN
use_release_issue()

# test
load_all()

adjust_uc(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "AGE",
  u_model_coefs = c(0.10, 0.10, 0.10, 0.02),
)

adjust_emc(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "AGE",
  x_model_coefs = c(0.10, 0.10, 0.10, 0.02),
)

set.seed(1)
adjust_sel(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "AGE",
  s_model_coefs = c(0.05, 0.20, 0.20)
)
# 2.28 (1.05, 4.96)