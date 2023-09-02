library(devtools)
library(tidyverse)
library(available)
library(roxygen2)

install.packages("languageserver")

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